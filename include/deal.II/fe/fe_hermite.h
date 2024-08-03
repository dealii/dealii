// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_hermite_h
#define dealii_fe_hermite_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_hermite.h>

#include <deal.II/fe/fe_poly.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup fe
 * @{
 */



/**
 * This class implements a Hermite interpolation basis of maximum regularity
 * elements (see @cite CiarletRiavart1972interpolation). These bases are always
 * of odd polynomial degree, have regularity $r=\frac{p-1}{2}$ and are defined
 * up to polynomial degree $p=13$, with larger degrees currently being
 * ill-conditioned.
 *
 * Each node has $(r+1)^{d}$ degrees of freedom (DoFs) assigned to it,
 * corresponding to various derivatives up to order $r$ in each direction.
 * DoFs at each node are not consecutive in lexicographic ordering for
 * $d>1$ due to the tensor product construction of the basis. The
 * ordering is determined by the direction of the derivative each function
 * corresponds to; first by $x$-derivatives, then $y$, then $z$. Locally over
 * each element the DoFs are ordered similarly. See below for the local ordering
 * for $r=1$, where DoFs are ordered from 0 to $(2r+2)^{d}-1$:
 *
 * <code>FE_Hermite<1>(3)</code>
 *
 * @verbatim
 * (0)________________(2)
 * (1)                (3)
 * @endverbatim
 *
 * <code>FE_Hermite<2>(3)</code>:
 *
 * @verbatim
 * ( 8, 9)__________(10,11)
 * (12,13)          (14,15)
 *    |                |
 *    |                |
 *    |                |
 *    |                |
 *    |                |
 *    |                |
 * ( 0, 1)__________( 2, 3)
 * ( 4, 5)          ( 6, 7)
 * @endverbatim
 *
 * <code>FE_Hermite<3>(3)</code>:
 *
 * @verbatim
 *        (40,41,44,45)__(42,43,46,47)          (40,41,44,45)__(42,43,46,47)
 *        (56,57,60,61)  (58,59,62,63)          (56,57,60,61)  (58,59,62,63)
 *           /  |              |                    /           /  |
 *          /   |              |                   /           /   |
 *         /    |              |                  /           /    |
 * (32,33,36,37)|              |         (32,33,36,37)__(34,35,38,39)
 * (48,49,52,53)|              |         (48,49,52,53)  (50,51,54,55)
 *        |     |              |                |            |     |
 *        |( 8,9,12,13 )__(10,11,14,15)         |            |(10,11,14,15)
 *        |(24,25,28,29)  (26,27,30,31)         |            |(26,27,30,31)
 *        |    /              /                 |            |    /
 *        |   /              /                  |            |   /
 *        |  /              /                   |            |  /
 *   (  0,1,4,5  )___(  2,3,6,7  )      (  0,1,4,5  )___(  2,3,6,7  )
 *   (16,17,20,21)   (18,19,22,23)      (16,17,20,21)   (18,19,22,23)
 * @endverbatim
 *
 * Due to the increased regularity constraints many of these functions are
 * shared between elements, leading to fewer DoFs overall.
 */
template <int dim, int spacedim = dim>
class FE_Hermite : public FE_Poly<dim, spacedim>
{
public:
  /**
   * Constructor that creates an Hermite finite element object with
   * continuity in derivatives up to and including order
   * <code>(fe_degree - 1)/2</code>.
   * Throws an exception if @p fe_degree is not an odd positive integer.
   */
  FE_Hermite(const unsigned int fe_degree);

  // Other functions
  /**
   * Returns <code>FE_Hermite<dim,spacedim>(fe_degree)</code> as a
   * <code>std::string</code> with @p dim, @p spacedim and @p fe_degree
   * replaced with the correct values.
   */
  virtual std::string
  get_name() const override;

  /**
   * @copydoc dealii::FiniteElement::clone()
   */
  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * @copydoc dealii::FiniteElement::requires_update_flags()
   */
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  /**
   * @copydoc dealii::FiniteElement::hp_vertex_dof_identities()
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * @copydoc dealii::FiniteElement::hp_line_dof_identities()
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * @copydoc dealii::FiniteElement::hp_quad_dof_identities()
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  /**
   * @copydoc FiniteElement::compare_for_domination()
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &other_fe,
                         const unsigned int codim) const override;

  /**
   * Returns the mapping between lexicographic and hierarchic numbering
   * schemes for Hermite. See the class documentation for diagrams of
   * examples of lexicographic numbering for Hermite elements.
   */
  std::vector<unsigned int>
  get_lexicographic_to_hierarchic_numbering() const;

  /**
   * This re-implements FiniteElement::fill_fe_values() for a Hermite
   * polynomial basis, to account for the more complicated shape function
   * re-scaling that a Hermite basis requires. On a non-uniform grid,
   * failing to rescale function values will lead to discontinuities across
   * element boundaries in derivatives if elements have different sizes.
   *
   * At present this is done by only allowing certain cell mappings
   * (currently only MappingCartesian) which guarantee
   * rectangular cells, which ensures that the directions of derivatives
   * does not change. The derivative continuity can then be enforced by
   * re-scaling shape functions so that the magnitude of each derivative is
   * equivalent to the derivative on the reference interval.
   */
  virtual void
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim>                                      &quadrature,
    const Mapping<dim, spacedim>                               &mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
                                                                  &mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  using FiniteElement<dim, spacedim>::fill_fe_face_values;

  /**
   * This re-implements FiniteElement::fill_fe_face_values() for
   * Hermite polynomial bases, for the same reasons as described for
   * FEHermite::fill_fe_values().
   */
  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1>                             &quadrature,
    const Mapping<dim, spacedim>                               &mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
                                                                  &mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;



protected:
  /**
   * Wrapper function for FE_Poly::get_data() that ensures relevant
   * shape value data is copied to the new data object as well.
   */
  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim>        &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override
  {
    std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
          data_ptr = FE_Poly<dim, spacedim>::get_data(update_flags,
                                                  mapping,
                                                  quadrature,
                                                  output_data);
    auto &data =
      dynamic_cast<typename FE_Poly<dim, spacedim>::InternalData &>(*data_ptr);

    const unsigned int n_q_points = quadrature.size();

    if ((update_flags & update_values) &&
        ((output_data.shape_values.n_rows() > 0) &&
         (output_data.shape_values.n_cols() == n_q_points)))
      data.shape_values = output_data.shape_values;

    return data_ptr;
  }



private:
  /**
   * Variable storing the order of the highest derivative that the current
   * FE_Hermite object can enforce continuity for. Here the order
   * of derivative only counts in one spatial direction, so the derivative
   * $\frac{d^{2}f}{dx \; dy}$ would be counted as a first order derivative
   * of $f$, as an example.
   */
  unsigned int regularity;



public:
  /**
   * Returns the regularity of the Hermite FE basis.
   */
  inline unsigned int
  get_regularity() const
  {
    return this->regularity;
  };
};

/** @} */



DEAL_II_NAMESPACE_CLOSE

#endif
