// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2014 by the deal.II authors
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

#ifndef __deal2__fe_face_h
#define __deal2__fe_face_h

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
 * Although the name does not give it away, the element is discontinuous
 * at locations where faces of cells meet. In particular,
 * this finite element is the trace space of FE_RaviartThomas on the faces and
 * serves in hybridized methods, e.g. in combination with the FE_DGQ
 * element. Its use is demonstrated in the step-51 tutorial program.
 *
 * @note Since this element is defined only on faces, only
 * FEFaceValues and FESubfaceValues will be able to extract reasonable
 * values from any face polynomial. In order to make the use of
 * FESystem simpler, using a (cell) FEValues object will not fail using this finite
 * element space, but all shape function values extracted will be equal
 * to zero.
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

  virtual FiniteElement<dim,spacedim> *clone() const;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_FaceQ<dim>(degree)</tt>, with <tt>dim</tt> and
   * <tt>degree</tt> replaced by appropriate values.
   */
  virtual std::string get_name () const;

  /**
   * Return the matrix interpolating from a face of of one element to the face
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
   * Return the matrix interpolating from a face of of one element to the face
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
   * For a definition of domination, see FiniteElementBase::Domination and in
   * particular the @ref hp_paper "hp paper".
   */
  virtual
  FiniteElementDomination::Domination
  compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * Returns a list of constant modes of the element. For this element, it
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

  /**
   * Clone method.
   */
  virtual FiniteElement<1,spacedim> *clone() const;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_FaceQ<dim>(degree)</tt>, with <tt>dim</tt> and
   * <tt>degree</tt> replaced by appropriate values.
   */
  virtual std::string get_name () const;

  /**
   * Return the matrix interpolating from a face of of one element to the face
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
   * Return the matrix interpolating from a face of of one element to the face
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
   * Return whether this element dominates the one given as argument when they
   * meet at a common face, whether it is the other way around, whether
   * neither dominates, or if either could dominate.
   *
   * For a definition of domination, see FiniteElementBase::Domination and in
   * particular the @ref hp_paper "hp paper".
   */
  virtual
  FiniteElementDomination::Domination
  compare_for_face_domination (const FiniteElement<1,spacedim> &fe_other) const;

  /**
   * Returns a list of constant modes of the element. For this element, it
   * simply returns one row with all entries set to true.
   */
  virtual std::pair<Table<2,bool>, std::vector<unsigned int> >
  get_constant_modes () const;

protected:
  virtual
  typename Mapping<1,spacedim>::InternalDataBase *
  get_data (const UpdateFlags,
            const Mapping<1,spacedim> &mapping,
            const Quadrature<1> &quadrature) const ;

  typename Mapping<1,spacedim>::InternalDataBase *
  get_face_data (const UpdateFlags,
                 const Mapping<1,spacedim> &mapping,
                 const Quadrature<0> &quadrature) const ;

  typename Mapping<1,spacedim>::InternalDataBase *
  get_subface_data (const UpdateFlags,
                    const Mapping<1,spacedim> &mapping,
                    const Quadrature<0> &quadrature) const ;

  virtual void
  fill_fe_values (const Mapping<1,spacedim>                           &mapping,
                  const typename Triangulation<1,spacedim>::cell_iterator &cell,
                  const Quadrature<1>                                 &quadrature,
                  typename Mapping<1,spacedim>::InternalDataBase      &mapping_internal,
                  typename Mapping<1,spacedim>::InternalDataBase      &fe_internal,
                  FEValuesData<1,spacedim>                            &data,
                  CellSimilarity::Similarity                       &cell_similarity) const;

  virtual void
  fill_fe_face_values (const Mapping<1,spacedim> &mapping,
                       const typename Triangulation<1,spacedim>::cell_iterator &cell,
                       const unsigned int                    face_no,
                       const Quadrature<0>                &quadrature,
                       typename Mapping<1,spacedim>::InternalDataBase      &mapping_internal,
                       typename Mapping<1,spacedim>::InternalDataBase      &fe_internal,
                       FEValuesData<1,spacedim> &data) const ;

  virtual void
  fill_fe_subface_values (const Mapping<1,spacedim> &mapping,
                          const typename Triangulation<1,spacedim>::cell_iterator &cell,
                          const unsigned int                    face_no,
                          const unsigned int                    sub_no,
                          const Quadrature<0>                &quadrature,
                          typename Mapping<1,spacedim>::InternalDataBase      &mapping_internal,
                          typename Mapping<1,spacedim>::InternalDataBase      &fe_internal,
                          FEValuesData<1,spacedim> &data) const ;


  /**
   * Determine the values that need to be computed on the unit cell to be able
   * to compute all values required by <tt>flags</tt>.
   *
   * For the purpuse of this function, refer to the documentation in
   * FiniteElement.
   *
   * This class assumes that shape functions of this FiniteElement do
   * <em>not</em> depend on the actual shape of the cells in real
   * space. Therefore, the effect in this element is as follows: if
   * <tt>update_values</tt> is set in <tt>flags</tt>, copy it to the
   * result. All other flags of the result are cleared, since everything else
   * must be computed for each cell.
   */
  virtual UpdateFlags update_once (const UpdateFlags flags) const;

  /**
   * Determine the values that need to be computed on every cell to be able to
   * compute all values required by <tt>flags</tt>.
   *
   * For the purpuse of this function, refer to the documentation in
   * FiniteElement.
   *
   * This class assumes that shape functions of this FiniteElement do
   * <em>not</em> depend on the actual shape of the cells in real space.
   *
   * The effect in this element is as follows:
   * <ul>
   *
   * <li> if <tt>update_gradients</tt> is set, the result will contain
   * <tt>update_gradients</tt> and <tt>update_covariant_transformation</tt>.
   * The latter is required to transform the gradient on the unit cell to the
   * real cell. Remark, that the action required by
   * <tt>update_covariant_transformation</tt> is actually performed by the
   * Mapping object used in conjunction with this finite element.
   *
   * <li> if <tt>update_hessians</tt> is set, the result will contain
   * <tt>update_hessians</tt> and <tt>update_covariant_transformation</tt>.
   * The rationale is the same as above and no higher derivatives of the
   * transformation are required, since we use difference quotients for the
   * actual computation.
   *
   * </ul>
   */
  virtual UpdateFlags update_each (const UpdateFlags flags) const;

private:
  /**
   * Return vector with dofs per vertex, line, quad, hex.
   */
  static std::vector<unsigned int> get_dpo_vector (const unsigned int deg);
};



/**
 * A finite element, which is a Legendre element of complete polynomials on
 * each face (i.e., it is the face equivalent of what FE_DGP is on cells)
 * and undefined in the interior of the cells. The basis functions on
 * the faces are from Polynomials::Legendre.
 *
 * Although the name does not give it away, the element is discontinuous
 * at locations where faces of cells meet. The element
 * serves in hybridized methods, e.g. in combination with the FE_DGP
 * element. An example of hybridizes methods can be found in the
 * step-51 tutorial program.
 *
 * @note Since this element is defined only on faces, only
 * FEFaceValues and FESubfaceValues will be able to extract reasonable
 * values from any face polynomial. In order to make the use of
 * FESystem simpler, using a (cell) FEValues object will not fail using this finite
 * element space, but all shape function values extracted will be equal
 * to zero.
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

  /**
   * Clone method.
   */
  virtual FiniteElement<dim,spacedim> *clone() const;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_FaceP<dim>(degree)</tt> , with <tt>dim</tt> and
   * <tt>degree</tt> replaced by appropriate values.
   */
  virtual std::string get_name () const;

  /**
   * Return the matrix interpolating from a face of of one element to the face
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
   * Return the matrix interpolating from a face of of one element to the face
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
   * For a definition of domination, see FiniteElementBase::Domination and in
   * particular the @ref hp_paper "hp paper".
   */
  virtual
  FiniteElementDomination::Domination
  compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * Returns a list of constant modes of the element. For this element, the
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
   * Returns the name of the element
   */
  std::string get_name() const;
};


DEAL_II_NAMESPACE_CLOSE

#endif
