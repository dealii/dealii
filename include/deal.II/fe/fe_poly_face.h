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

#ifndef dealii_fe_poly_face_h
#define dealii_fe_poly_face_h


#include <deal.II/base/config.h>

#include <deal.II/base/qprojector.h>

#include <deal.II/fe/fe.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup fe
 * @{
 */

/**
 * @warning This class is not sufficiently tested yet!
 *
 * This class gives a unified framework for the implementation of
 * FiniteElement classes only located on faces of the mesh. They are based on
 * polynomial spaces like the TensorProductPolynomials or a PolynomialSpace
 * classes.
 *
 * Every class that implements the following functions can be used as template
 * parameter PolynomialType.
 *
 * @code
 * double compute_value (const unsigned int i,
 *                       const Point<dim> &p) const;
 * @endcode
 * Example classes are TensorProductPolynomials, PolynomialSpace or
 * PolynomialsP.
 *
 * This class is not a fully implemented FiniteElement class. Instead there
 * are several pure virtual functions declared in the FiniteElement class
 * which cannot be implemented by this class but are left for implementation
 * in derived classes.
 */
template <typename PolynomialType,
          int dim      = PolynomialType::dimension + 1,
          int spacedim = dim>
class FE_PolyFace : public FiniteElement<dim, spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_PolyFace(const PolynomialType         &poly_space,
              const FiniteElementData<dim> &fe_data,
              const std::vector<bool>      &restriction_is_additive_flags);

  /**
   * Return the polynomial degree of this finite element, i.e. the value
   * passed to the constructor.
   */
  unsigned int
  get_degree() const;

  // for documentation, see the FiniteElement base class
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

protected:
  /*
   * NOTE: The following functions have their definitions inlined into the class
   * declaration because we otherwise run into a compiler error with MS Visual
   * Studio.
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
    (void)update_flags;
    (void)mapping;
    (void)quadrature;
    (void)output_data;
    return std::make_unique<InternalData>();
  }

  using FiniteElement<dim, spacedim>::get_face_data;

  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_face_data(
    const UpdateFlags               update_flags,
    const Mapping<dim, spacedim>   &mapping,
    const hp::QCollection<dim - 1> &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override
  {
    (void)mapping;
    (void)output_data;
    AssertDimension(quadrature.size(), 1);

    // generate a new data object and
    // initialize some fields
    std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
          data_ptr   = std::make_unique<InternalData>();
    auto &data       = dynamic_cast<InternalData &>(*data_ptr);
    data.update_each = requires_update_flags(update_flags);

    const unsigned int n_q_points = quadrature[0].size();

    // some scratch arrays
    std::vector<double>             values(0);
    std::vector<Tensor<1, dim - 1>> grads(0);
    std::vector<Tensor<2, dim - 1>> grad_grads(0);
    std::vector<Tensor<3, dim - 1>> empty_vector_of_3rd_order_tensors;
    std::vector<Tensor<4, dim - 1>> empty_vector_of_4th_order_tensors;

    // initialize fields only if really
    // necessary. otherwise, don't
    // allocate memory
    if (data.update_each & update_values)
      {
        values.resize(poly_space.n());
        data.shape_values.resize(poly_space.n(),
                                 std::vector<double>(n_q_points));
        for (unsigned int i = 0; i < n_q_points; ++i)
          {
            poly_space.evaluate(quadrature[0].point(i),
                                values,
                                grads,
                                grad_grads,
                                empty_vector_of_3rd_order_tensors,
                                empty_vector_of_4th_order_tensors);

            for (unsigned int k = 0; k < poly_space.n(); ++k)
              data.shape_values[k][i] = values[k];
          }
      }
    // No derivatives of this element
    // are implemented.
    if (data.update_each & update_gradients ||
        data.update_each & update_hessians)
      {
        DEAL_II_NOT_IMPLEMENTED();
      }

    return data_ptr;
  }

  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_subface_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim - 1>    &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override
  {
    return get_face_data(
      update_flags,
      mapping,
      hp::QCollection<dim - 1>(QProjector<dim - 1>::project_to_all_children(
        ReferenceCells::get_hypercube<dim - 1>(), quadrature)),
      output_data);
  }

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
   * Fields of cell-independent data.
   *
   * For information about the general purpose of this class, see the
   * documentation of the base class.
   */
  class InternalData : public FiniteElement<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * Array with shape function values in quadrature points on one face.
     * There is one row for each shape function, containing values for each
     * quadrature point.
     *
     * In this array, we store the values of the shape function in the
     * quadrature points on one face of the unit cell. Since these values do
     * not change under transformation to the real cell, we only need to copy
     * them over when visiting a concrete cell.
     *
     * In particular, we can simply copy the same set of values to each of the
     * faces.
     */
    std::vector<std::vector<double>> shape_values;
  };

  /**
   * The polynomial space. Its type is given by the template parameter
   * PolynomialType.
   */
  PolynomialType poly_space;
};

/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif
