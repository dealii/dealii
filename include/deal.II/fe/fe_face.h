// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2017 by the deal.II authors
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

#ifndef dealii_fe_face_h
#define dealii_fe_face_h

#include <deal.II/base/config.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/fe/fe_poly_face.h>

DEAL_II_NAMESPACE_OPEN


/**
 * A finite element, which is a tensor product polynomial on each face and
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
 * FESubfaceValues will be able to extract reasonable values from any face
 * polynomial. In order to make the use of FESystem simpler, using a (cell)
 * FEValues object will not fail using this finite element space, but all
 * shape function values extracted will be equal to zero.
 *
 * @ingroup fe
 * @author Guido Kanschat, Martin Kronbichler
 * @date 2009, 2011, 2013
 */
template <int dim, int spacedim=dim>
class FE_FaceQ : public FE_PolyFace<TensorProductPolynomials<dim-1>, dim, spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials of degree <tt>p</tt>. The
   * shape functions created using this constructor correspond to Lagrange
   * polynomials in each coordinate direction.
   */
  FE_FaceQ (const unsigned int p);

  virtual
  std::unique_ptr<FiniteElement<dim,spacedim> >
  clone() const;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_FaceQ<dim>(degree)</tt>, with <tt>dim</tt> and
   * <tt>degree</tt> replaced by appropriate values.
   */
  virtual std::string get_name () const;

  /**
   * Implementation of the corresponding function in the FiniteElement
   * class.  Since the current element is interpolatory, the nodal
   * values are exactly the support point values. Furthermore, since
   * the current element is scalar, the support point values need to
   * be vectors of length 1.
   */
  virtual
  void
  convert_generalized_support_point_values_to_dof_values (const std::vector<Vector<double> > &support_point_values,
                                                          std::vector<double>                &nodal_values) const;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                 FullMatrix<double>       &matrix) const;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                    const unsigned int        subface,
                                    FullMatrix<double>       &matrix) const;

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

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
   * The set of such constraints is non-empty only for dim==1.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_vertex_dof_identities (const FiniteElement<dim, spacedim> &fe_other) const;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   *
   * The set of such constraints is non-empty only for dim==2.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_line_dof_identities (const FiniteElement<dim, spacedim> &fe_other) const;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on quads.
   *
   * The set of such constraints is non-empty only for dim==3.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_quad_dof_identities (const FiniteElement<dim, spacedim> &fe_other) const;

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp compatible".
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
   * Return a list of constant modes of the element. For this element, it
   * simply returns one row with all entries set to true.
   */
  virtual std::pair<Table<2,bool>, std::vector<unsigned int> >
  get_constant_modes () const;

private:
  /**
   * Return vector with dofs per vertex, line, quad, hex.
   */
  static std::vector<unsigned int> get_dpo_vector (const unsigned int deg);
};



/**
 * Specialization of FE_FaceQ for 1D. In that case, the finite element only
 * consists of one degree of freedom in each of the two faces (= vertices) of
 * a cell, irrespective of the degree. However, this element still accepts a
 * degree in its constructor and also returns that degree. This way,
 * dimension-independent programming with trace elements is also possible in
 * 1D (even though there is no computational benefit at all from it in 1D).
 *
 * @ingroup fe
 * @author Guido Kanschat, Martin Kronbichler
 * @date 2014
 */
template <int spacedim>
class FE_FaceQ<1,spacedim> : public FiniteElement<1,spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_FaceQ (const unsigned int p);

  virtual
  std::unique_ptr<FiniteElement<1,spacedim>>
                                          clone() const;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_FaceQ<dim>(degree)</tt>, with <tt>dim</tt> and
   * <tt>degree</tt> replaced by appropriate values.
   */
  virtual std::string get_name () const;

  // for documentation, see the FiniteElement base class
  virtual
  UpdateFlags
  requires_update_flags (const UpdateFlags update_flags) const;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_face_interpolation_matrix (const FiniteElement<1,spacedim> &source,
                                 FullMatrix<double>       &matrix) const;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_subface_interpolation_matrix (const FiniteElement<1,spacedim> &source,
                                    const unsigned int        subface,
                                    FullMatrix<double>       &matrix) const;

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp compatible".
   */
  virtual bool hp_constraints_are_implemented () const;

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
   * The set of such constraints is non-empty only for dim==1.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_vertex_dof_identities (const FiniteElement<1, spacedim> &fe_other) const;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   *
   * The set of such constraints is non-empty only for dim==2.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_line_dof_identities (const FiniteElement<1, spacedim> &fe_other) const;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on quads.
   *
   * The set of such constraints is non-empty only for dim==3.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_quad_dof_identities (const FiniteElement<1, spacedim> &fe_other) const;

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
  compare_for_face_domination (const FiniteElement<1,spacedim> &fe_other) const;

  /**
   * Return a list of constant modes of the element. For this element, it
   * simply returns one row with all entries set to true.
   */
  virtual std::pair<Table<2,bool>, std::vector<unsigned int> >
  get_constant_modes () const;

protected:
  /*
   * NOTE: The following functions have their definitions inlined into the class declaration
   * because we otherwise run into a compiler error with MS Visual Studio.
   */


  virtual
  std::unique_ptr<typename FiniteElement<1,spacedim>::InternalDataBase>
  get_data (const UpdateFlags                                                  /*update_flags*/,
            const Mapping<1,spacedim>                                         &/*mapping*/,
            const Quadrature<1>                                               &/*quadrature*/,
            dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1, spacedim> &/*output_data*/) const
  {
    return std_cxx14::make_unique<typename FiniteElement<1, spacedim>::InternalDataBase>();
  }

  std::unique_ptr<typename FiniteElement<1,spacedim>::InternalDataBase>
  get_face_data(const UpdateFlags update_flags,
                const Mapping<1,spacedim> &/*mapping*/,
                const Quadrature<0> &quadrature,
                dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1, spacedim> &/*output_data*/) const
  {
    // generate a new data object and initialize some fields
    auto data = std_cxx14::make_unique<typename FiniteElement<1,spacedim>::InternalDataBase>();
    data->update_each = requires_update_flags(update_flags);

    const unsigned int n_q_points = quadrature.size();
    AssertDimension(n_q_points, 1);
    (void)n_q_points;

    // No derivatives of this element are implemented.
    if (data->update_each & update_gradients || data->update_each & update_hessians)
      {
        Assert(false, ExcNotImplemented());
      }

    return std::move(data);
  }

  std::unique_ptr<typename FiniteElement<1,spacedim>::InternalDataBase>
  get_subface_data(const UpdateFlags                                                  update_flags,
                   const Mapping<1,spacedim>                                         &mapping,
                   const Quadrature<0>                                               &quadrature,
                   dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1, spacedim> &output_data) const
  {
    return get_face_data(update_flags, mapping, quadrature, output_data);
  }

  virtual
  void
  fill_fe_values (const typename Triangulation<1,spacedim>::cell_iterator           &cell,
                  const CellSimilarity::Similarity                                   cell_similarity,
                  const Quadrature<1>                                               &quadrature,
                  const Mapping<1,spacedim>                                         &mapping,
                  const typename Mapping<1,spacedim>::InternalDataBase              &mapping_internal,
                  const dealii::internal::FEValuesImplementation::MappingRelatedData<1, spacedim> &mapping_data,
                  const typename FiniteElement<1,spacedim>::InternalDataBase        &fe_internal,
                  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1, spacedim> &output_data) const;

  virtual
  void
  fill_fe_face_values (const typename Triangulation<1,spacedim>::cell_iterator           &cell,
                       const unsigned int                                                 face_no,
                       const Quadrature<0>                                               &quadrature,
                       const Mapping<1,spacedim>                                         &mapping,
                       const typename Mapping<1,spacedim>::InternalDataBase              &mapping_internal,
                       const dealii::internal::FEValuesImplementation::MappingRelatedData<1, spacedim> &mapping_data,
                       const typename FiniteElement<1,spacedim>::InternalDataBase        &fe_internal,
                       dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1, spacedim> &output_data) const;

  virtual
  void
  fill_fe_subface_values (const typename Triangulation<1,spacedim>::cell_iterator           &cell,
                          const unsigned int                                                 face_no,
                          const unsigned int                                                 sub_no,
                          const Quadrature<0>                                               &quadrature,
                          const Mapping<1,spacedim>                                         &mapping,
                          const typename Mapping<1,spacedim>::InternalDataBase              &mapping_internal,
                          const dealii::internal::FEValuesImplementation::MappingRelatedData<1, spacedim> &mapping_data,
                          const typename FiniteElement<1,spacedim>::InternalDataBase        &fe_internal,
                          dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1, spacedim> &output_data) const;

private:
  /**
   * Return vector with dofs per vertex, line, quad, hex.
   */
  static
  std::vector<unsigned int>
  get_dpo_vector (const unsigned int deg);
};



/**
 * A finite element, which is a Legendre element of complete polynomials on
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
 * FESubfaceValues will be able to extract reasonable values from any face
 * polynomial. In order to make the use of FESystem simpler, using a (cell)
 * FEValues object will not fail using this finite element space, but all
 * shape function values extracted will be equal to zero.
 *
 * @ingroup fe
 * @author Martin Kronbichler
 * @date 2013
 */
template <int dim, int spacedim=dim>
class FE_FaceP : public FE_PolyFace<PolynomialSpace<dim-1>, dim, spacedim>
{
public:
  /**
   * Constructor for complete basis of polynomials of degree <tt>p</tt>. The
   * shape functions created using this constructor correspond to Legendre
   * polynomials in each coordinate direction.
   */
  FE_FaceP(unsigned int p);

  virtual
  std::unique_ptr<FiniteElement<dim,spacedim> >
  clone() const;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_FaceP<dim>(degree)</tt> , with <tt>dim</tt> and
   * <tt>degree</tt> replaced by appropriate values.
   */
  virtual std::string get_name () const;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                 FullMatrix<double>       &matrix) const;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. This
   * element only provides interpolation matrices for elements of the same
   * type and FE_Nothing. For all other elements, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                    const unsigned int        subface,
                                    FullMatrix<double>       &matrix) const;

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp compatible".
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
   * Return a list of constant modes of the element. For this element, the
   * first entry on each face is true, all other are false (as the constant
   * function is represented by the first base function of Legendre
   * polynomials).
   */
  virtual std::pair<Table<2,bool>, std::vector<unsigned int> >
  get_constant_modes () const;

private:
  /**
   * Return vector with dofs per vertex, line, quad, hex.
   */
  static std::vector<unsigned int> get_dpo_vector (const unsigned int deg);
};



/**
 * FE_FaceP in 1D, i.e., with degrees of freedom on the element vertices.
 */
template <int spacedim>
class FE_FaceP<1,spacedim> : public FE_FaceQ<1,spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_FaceP (const unsigned int p);

  /**
   * Return the name of the element
   */
  std::string get_name() const;
};


DEAL_II_NAMESPACE_CLOSE

#endif
