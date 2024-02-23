// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_nothing_h
#define dealii_fe_nothing_h

#include <deal.II/base/config.h>

#include <deal.II/fe/fe.h>

DEAL_II_NAMESPACE_OPEN


/**
 * @addtogroup fe
 * @{
 */

/**
 * Definition of a finite element space with zero degrees of freedom and that,
 * consequently, can only represent a single function: the zero function.
 *
 * This class is useful (in the context of an hp-method) to represent empty
 * cells in the triangulation on which no degrees of freedom should be
 * allocated, or to describe a field that is extended by zero to a part of the
 * domain where we don't need it. Thus a triangulation may be divided into two
 * regions: an active region where normal elements are used, and an inactive
 * region where FE_Nothing elements are used. The DoFHandler will therefore
 * assign no degrees of freedom to the FE_Nothing cells, and this subregion is
 * therefore implicitly deleted from the computation. step-10 and step-46 show
 * use cases for this element. An interesting application for this element is
 * also presented in the paper @cite Cangiani2012.
 *
 *
 * <h3>FE_Nothing as seen as a function space</h3>
 *
 * Finite elements are often best interpreted as forming a
 * [function space](https://en.wikipedia.org/wiki/Function_space), i.e., a
 * set of functions that form a
 * [vector space](https://en.wikipedia.org/wiki/Vector_space). One can indeed
 * interpret FE_Nothing in this light: It corresponds to the function space
 * $V_h=\{0\}$, i.e., the set of functions that are zero everywhere.
 * (The constructor can take an argument that, if greater than one, extends
 * the space to one of vector-valued functions with more than one component,
 * with all components equal to zero everywhere.) Indeed, this is a vector
 * space since every linear combination of elements in the vector space is
 * also an element in the vector space, as is every multiple of the single
 * element zero. It is obvious that the function space has no degrees of
 * freedom, thus the name of the class.
 *
 *
 * <h3>FE_Nothing in combination with other elements</h3>
 *
 * In situations such as those of step-46, one uses FE_Nothing on cells
 * where one is not interested in a solution variable. For example, in fluid
 * structure interaction problems, the fluid velocity is only defined on
 * cells inside the fluid part of the domain. One then uses FE_Nothing
 * on cells in the solid part of the domain to describe the finite element
 * space for the velocity. In other words, the velocity lives everywhere
 * conceptually, but it is identically zero in those parts of the domain
 * where it is not of interest and doesn't use up any degrees of freedom
 * there.
 *
 * The question is what happens at the interface between areas where one
 * is interested in the solution (and uses a "normal" finite element) and
 * where one is not interested (and uses FE_Nothing): Should the solution
 * at that interface be zero -- i.e., we consider a "continuous" finite
 * element field that happens to be zero in that area where FE_Nothing
 * is used -- or is there no requirement for continuity at the interface.
 * In the deal.II language, this is encoded by what the function
 * FiniteElement::compare_for_domination() returns: If the FE_Nothing
 * "dominates", then the solution must be zero at the interface; if it
 * does not, then there is no requirement and one can think of FE_Nothing
 * as a function space that is in general discontinuous (i.e., there is
 * no requirement for any kind of continuity at cell interfaces) but on
 * every cell equal to zero.
 *
 * A constructor argument denotes whether the element should be considered
 * dominating or not. The default is for it not to dominate, i.e.,
 * FE_Nothing is treated as a discontinuous element.
 *
 *
 * <h3>FE_Nothing in the context of hanging nodes</h3>
 *
 * Note that some care must be taken that the resulting mesh topology
 * continues to make sense when FE_Nothing elements are introduced. This is
 * particularly true when dealing with hanging node constraints, because the
 * library makes some basic assumptions about the nature of those constraints.
 * The following geometries are acceptable:
 * @code
 * +---------+----+----+
 * |         | 0  |    |
 * |    1    +----+----+
 * |         | 0  |    |
 * +---------+----+----+
 * @endcode
 * @code
 * +---------+----+----+
 * |         | 1  |    |
 * |    0    +----+----+
 * |         | 1  |    |
 * +---------+----+----+
 * @endcode
 * Here, 0 denotes an FE_Nothing cell, and 1 denotes some other element type.
 * The library has no difficulty computing the necessary hanging node
 * constraints in these cases (i.e. no constraint). However, the following
 * geometry is NOT acceptable (at least in the current implementation):
 * @code
 * +---------+----+----+
 * |         | 0  |    |
 * |    1    +----+----+
 * |         | 1  |    |
 * +---------+----+----+
 * @endcode
 * The distinction lies in the mixed nature of the child faces, a case we have
 * not implemented as of yet.
 */
template <int dim, int spacedim = dim>
class FE_Nothing : public FiniteElement<dim, spacedim>
{
public:
  /**
   * Constructor.
   *
   * @param[in] type Specifies the reference-cell type.
   *
   * @param[in] n_components Denotes the number of
   * vector components to give this finite element. The default is one.
   *
   * @param[in] dominate Decides whether FE_Nothing will dominate
   * any other FE in compare_for_domination() (with the default being `false`).
   * Therefore at interfaces where, for example, a $Q_1$ meets an FE_Nothing, we
   * will force the traces of the two functions to be the same. Because the
   * FE_Nothing encodes a space that is zero everywhere, this means that the
   * $Q_1$ field will be forced to become zero at this interface. See also the
   * discussion in the general documentation of this class.
   */
  FE_Nothing(const ReferenceCell &type,
             const unsigned int   n_components = 1,
             const bool           dominate     = false);

  /**
   * Same as above but for a hypercube reference-cell type.
   */
  FE_Nothing(const unsigned int n_components = 1, const bool dominate = false);

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * Return a string that uniquely identifies a finite element. The name is
   * <tt>FE_Nothing@<dim,spacedim@>(type, n_components, dominating)</tt> where
   * <tt>dim</tt>, <tt>spacedim</tt>, <tt>type</tt>, and <tt>n_components</tt>
   * are all specified by the constructor or type signature with the following
   * exceptions:
   * <ol>
   *   <li>If <tt>spacedim == dim</tt> then that field is not printed.</li>
   *   <li>If <tt>type</tt> is a hypercube then that field is not printed.</li>
   *   <li>If <tt>n_components == 1</tt> then that field is not printed.</li>
   *   <li>If <tt>dominate == false</tt> then that field is not printed.</li>
   * </ol>
   */
  virtual std::string
  get_name() const override;

  // for documentation, see the FiniteElement base class
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  /**
   * Return the value of the @p ith shape function at the point @p p. @p p is
   * a point on the reference element. Because the current element has no
   * degrees of freedom, this function should obviously not be called in
   * practice.  All this function really does, therefore, is trigger an
   * exception.
   */
  virtual double
  shape_value(const unsigned int i, const Point<dim> &p) const override;

  virtual void
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim>                                      &quadrature,
    const Mapping<dim, spacedim>                               &mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                                  &mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  using FiniteElement<dim, spacedim>::fill_fe_face_values;

  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1>                             &quadrature,
    const Mapping<dim, spacedim>                               &mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                                  &mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          sub_no,
    const Quadrature<dim - 1>                                  &quadrature,
    const Mapping<dim, spacedim>                               &mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                                  &mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  /**
   * Prepare internal data structures and fill in values independent of the
   * cell. Returns a pointer to an object of which the caller of this function
   * then has to assume ownership (which includes destruction when it is no
   * more needed).
   *
   * In the current case, this function just returns a default pointer, since
   * no meaningful data exists for this element.
   */
  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim>        &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  /**
   * @copydoc FiniteElement::compare_for_domination()
   *
   * In the current case, this element is assumed to dominate if the second
   * argument in the constructor @p dominate is true. When this argument is
   * false and @p fe_other is also of type FE_Nothing(), either element can
   * dominate. Otherwise there are no_requirements.
   *
   * See also the discussion in the general documentation of this class.
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;



  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * Return the matrix interpolating from the given finite element to the
   * present one. Since the current finite element has no degrees of freedom,
   * the interpolation matrix is necessarily empty.
   */
  virtual void
  get_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source_fe,
    FullMatrix<double>                 &interpolation_matrix) const override;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element. The size of the matrix is then
   * <tt>source.#dofs_per_face</tt> times <tt>this->#dofs_per_face</tt>.
   *
   * Since the current finite element has no degrees of freedom, the
   * interpolation matrix is necessarily empty.
   */

  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source_fe,
                                FullMatrix<double> &interpolation_matrix,
                                const unsigned int  face_no = 0) const override;


  /**
   * Return the matrix interpolating from a face of one element to the
   * subface of the neighboring element. The size of the matrix is then
   * <tt>source.#dofs_per_face</tt> times <tt>this->#dofs_per_face</tt>.
   *
   * Since the current finite element has no degrees of freedom, the
   * interpolation matrix is necessarily empty.
   */

  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source_fe,
    const unsigned int                  index,
    FullMatrix<double>                 &interpolation_matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * Return a list of constant modes of the element.
   *
   * Since the current finite element has no degrees of freedom, the returned
   * list is necessarily empty.
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  /**
   * @return true if the FE dominates any other.
   */
  bool
  is_dominating() const;

private:
  /**
   * If true, this element will dominate any other apart from itself in
   * compare_for_domination(). This is because a space that only contains the
   * zero function is definitely smaller (and consequently dominant) when
   * compared to any other finite element space.
   */
  const bool dominate;
};


/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif
