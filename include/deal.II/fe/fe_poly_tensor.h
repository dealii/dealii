// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_poly_tensor_h
#define dealii_fe_poly_tensor_h


#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/mutex.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/tensor_polynomials_base.h>

#include <deal.II/fe/fe.h>

#include <deal.II/lac/full_matrix.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

/**
 * This class provides a unified framework for the implementation of
 * FiniteElement classes based on tensor-valued polynomial spaces like
 * PolynomialsBDM and PolynomialsRaviartThomas. In this, it is the
 * tensor-valued equivalent of the FE_Poly class.
 *
 * In essence, what this class requires is that a derived class describes
 * to it a (vector-valued) polynomial space in which every polynomial
 * has exactly @p dim vector components. The classes that provide such
 * implementations are all derived from the TensorPolynomialsBase
 * class, and an object of one of these derived types needs to be
 * provided to the constructor of this class.
 *
 *
 * <h3>Deriving classes</h3>
 *
 * This class is not a fully implemented FiniteElement class, but implements
 * some common features of vector valued elements based on vector valued
 * polynomial classes. What's missing here in particular is information on the
 * topological location of the node values (i.e., whether a degree of freedom
 * is logically at a vertex, edge, face, or the interior of a cell --
 * information that determines the continuity properties of the associated
 * shape functions across cell interfaces), and derived classes need to provide
 * this information.
 *
 * Similarly, in many cases, node functionals depend on the shape of the mesh
 * cell, since they evaluate normal or tangential components on the faces. In
 * order to allow for a set of transformations, the variable #mapping_kind has
 * been introduced. It needs be set in the constructor of a derived
 * class.
 *
 * Any derived class must decide on the polynomial space to use.  This
 * polynomial space should be implemented simply as a set of vector valued
 * polynomials like PolynomialsBDM and PolynomialsRaviartThomas.  In order to
 * facilitate this implementation, which basis the polynomial space chooses
 * is not of importance to the current class -- as described next, this class
 * handles the transformation from the basis chosen by the polynomial space
 * template argument to the basis we want to use for finite element
 * computations internally.
 *
 *
 * <h4>Determining the correct basis</h4>
 *
 * In most cases, the basis used by the class that describes the polynomial
 * space, $\{\tilde\varphi_j(\hat{\mathbf x})\}$, does not match the one we
 * want to use for the finite element description,
 * $\{\varphi_j(\hat{\mathbf x})\}$. Rather, we need to express the finite
 * element shape functions as a linear combination of the basis provided
 * by the polynomial space:
 * @f{align*}{
 *   \varphi_j = \sum_k c_{jk} \tilde\varphi_j.
 * @f}
 * These expansion coefficients $c_{jk}$ are typically computed in the
 * constructors of derived classes. To facilitate this, this class
 * at first (unless told otherwise, see below), assumes that the shape
 * functions should be exactly the ones provided by the polynomial
 * space. In the constructor of the derived class, one then typically has
 * code of the form
 * @code
 *   // Now compute the inverse node matrix, generating the correct
 *   // basis functions from the raw ones. For a discussion of what
 *   // exactly happens here, see FETools::compute_node_matrix.
 *   const FullMatrix<double> M = FETools::compute_node_matrix(*this);
 *   this->inverse_node_matrix.reinit(n_dofs, n_dofs);
 *   this->inverse_node_matrix.invert(M);
 *   // From now on, the shape functions provided by FiniteElement::shape_value
 *   // and similar functions will be the correct ones, not
 *   // the raw shape functions from the polynomial space anymore.
 * @endcode
 * The FETools::compute_node_matrix() function explains in more
 * detail what exactly it computes, and how; in any case, the result
 * is that @p inverse_node_matrix now contains the expansion coefficients
 * $c_{jk}$, and the fact that this block of code now sets the
 * matrix to a non-zero size indicates to the functions of the current
 * class that it should from then on use the expanded basis,
 * $\{\varphi_j(\hat{\mathbf x})\}$, and no longer the original, "raw"
 * basis $\{\tilde\varphi_j(\hat{\mathbf x})\}$ when asked for values
 * or derivatives of shape functions.
 *
 * In order for this scheme to work, it is important to ensure that
 * the size of the @p inverse_node_matrix be zero at the time when
 * FETools::compute_node_matrix() is called; thus, the call to this
 * function cannot be inlined into the last line -- the result of
 * the call really does need to be stored in the temporary object
 * @p M.
 *
 *
 * <h4>Setting the transformation</h4>
 *
 * In most cases, vector valued basis functions must be transformed when
 * mapped from the reference cell to the actual grid cell. These
 * transformations can be selected from the set MappingKind and stored in
 * #mapping_kind. Therefore, each constructor should contain a line like:
 * @code
 * this->mapping_kind = {mapping_none};
 * @endcode
 * (in case no mapping is required) or using whatever value among
 * the ones defined in MappingKind is appropriate for the element you
 * are implementing. If each shape function may be mapped by different
 * mappings, then @p mapping_kind may be a vector with the same number
 * of elements as there are shape functions.
 *
 * @see TensorPolynomialsBase
 * @ingroup febase
 */
template <int dim, int spacedim = dim>
class FE_PolyTensor : public FiniteElement<dim, spacedim>
{
public:
  /**
   * Constructor. This constructor does a deep copy of the polynomials
   * object via the TensorPolynomialsBase::clone() function and stores
   * a pointer to the copy. As a consequence, the calling site can
   * simply pass a temporary object as the first argument.
   */
  FE_PolyTensor(const TensorPolynomialsBase<dim> &polynomials,
                const FiniteElementData<dim>     &fe_data,
                const std::vector<bool>          &restriction_is_additive_flags,
                const std::vector<ComponentMask> &nonzero_components);


  /**
   * Copy constructor.
   */
  FE_PolyTensor(const FE_PolyTensor &fe);

  // for documentation, see the FiniteElement base class
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  /**
   * Compute the (scalar) value of shape function @p i at the given quadrature
   * point @p p. Since the elements represented by this class are vector
   * valued, there is no such scalar value and the function therefore throws
   * an exception.
   */
  virtual double
  shape_value(const unsigned int i, const Point<dim> &p) const override;

  // documentation inherited from the base class
  virtual double
  shape_value_component(const unsigned int i,
                        const Point<dim>  &p,
                        const unsigned int component) const override;

  /**
   * Compute the gradient of (scalar) shape function @p i at the given
   * quadrature point @p p. Since the elements represented by this class are
   * vector valued, there is no such scalar value and the function therefore
   * throws an exception.
   */
  virtual Tensor<1, dim>
  shape_grad(const unsigned int i, const Point<dim> &p) const override;

  // documentation inherited from the base class
  virtual Tensor<1, dim>
  shape_grad_component(const unsigned int i,
                       const Point<dim>  &p,
                       const unsigned int component) const override;

  /**
   * Compute the Hessian of (scalar) shape function @p i at the given
   * quadrature point @p p. Since the elements represented by this class are
   * vector valued, there is no such scalar value and the function therefore
   * throws an exception.
   */
  virtual Tensor<2, dim>
  shape_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  // documentation inherited from the base class
  virtual Tensor<2, dim>
  shape_grad_grad_component(const unsigned int i,
                            const Point<dim>  &p,
                            const unsigned int component) const override;

protected:
#ifndef DOXYGEN
  class InternalData;
#endif

  /**
   * Internal function called by fill_fe_values(), fill_fe_face_values(),
   * and fill_fe_subface_values().
   */
  void
  compute_fill(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const typename QProjector<dim>::DataSetDescriptor          &offset,
    const unsigned int                                          n_q_points,
    const Mapping<dim, spacedim>                               &mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                              &mapping_data,
    const typename FE_PolyTensor<dim, spacedim>::InternalData &fe_internal,
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
      &output_data) const;

  /**
   * The mapping type to be used to map shape functions from the reference
   * cell to the mesh cell. If this vector is length one, the same mapping
   * will be applied to all shape functions. If the vector size is equal to
   * the finite element dofs per cell, then each shape function will be mapped
   * according to the corresponding entry in the vector.
   */
  std::vector<MappingKind> mapping_kind;

  /**
   * Returns a boolean that is true when the finite element uses a single
   * mapping and false when the finite element uses multiple mappings.
   */
  bool
  single_mapping_kind() const;

  /**
   * For faces with non-standard face_orientation in 3d, the dofs on faces
   * (quads) have to be permuted in order to be combined with the correct shape
   * functions and additionally can change the sign. Given a local dof @p index
   * on a quad, return the sign of the permuted shape function.
   *
   * The permutation itself is returned by
   * adjust_quad_dof_index_for_face_orientation implemented in the interface
   * class FiniteElement<dim>.
   */
  bool
  adjust_quad_dof_sign_for_face_orientation(
    const unsigned int                 index,
    const unsigned int                 face_no,
    const types::geometric_orientation combined_orientation) const;

  /**
   * For faces with non-standard face_orientation in 3d, the dofs on faces
   * (quads) need not only to be permuted in order to be combined with the
   * correct shape functions. Additionally they may change their sign.
   *
   * The constructor of this class fills this table with 'false' values, i.e.,
   * no sign change at all. Derived finite element classes have to
   * fill this Table with the correct values, see the documentation in
   * GeometryInfo<dim> and this
   * @ref GlossCombinedOrientation "glossary entry on face orientations".
   *
   * The table must be filled in finite element classes derived
   * from FE_PolyTensor in a meaningful way since the permutation
   * pattern and the pattern of sign changes depends on how the finite element
   * distributes the local dofs on the faces. An example is the function
   * `initialize_quad_dof_index_permutation_and_sign_change()` in the
   * FE_RaviartThomas class that fills this table.
   */
  std::vector<Table<2, bool>> adjust_quad_dof_sign_for_face_orientation_table;

  /**
   * Returns MappingKind @p i for the finite element.
   */
  MappingKind
  get_mapping_kind(const unsigned int i) const;

  /* NOTE: The following function has its definition inlined into the class
     declaration because we otherwise run into a compiler error with MS Visual
     Studio. */
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
    (void)mapping;
    (void)output_data;
    // generate a new data object and
    // initialize some fields
    std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
          data_ptr   = std::make_unique<InternalData>();
    auto &data       = dynamic_cast<InternalData &>(*data_ptr);
    data.update_each = requires_update_flags(update_flags);

    const unsigned int n_q_points = quadrature.size();

    // some scratch arrays
    std::vector<Tensor<1, dim>> values(0);
    std::vector<Tensor<2, dim>> grads(0);
    std::vector<Tensor<3, dim>> grad_grads(0);
    std::vector<Tensor<4, dim>> third_derivatives(0);
    std::vector<Tensor<5, dim>> fourth_derivatives(0);

    if (update_flags & (update_values | update_gradients | update_hessians))
      data.dof_sign_change.resize(this->dofs_per_cell);

    // initialize fields only if really
    // necessary. otherwise, don't
    // allocate memory

    const bool update_transformed_shape_values =
      std::any_of(this->mapping_kind.begin(),
                  this->mapping_kind.end(),
                  [](const MappingKind t) { return t != mapping_none; });

    const bool update_transformed_shape_grads =
      std::any_of(this->mapping_kind.begin(),
                  this->mapping_kind.end(),
                  [](const MappingKind t) {
                    return (t == mapping_raviart_thomas || t == mapping_piola ||
                            t == mapping_nedelec || t == mapping_contravariant);
                  });

    const bool update_transformed_shape_hessian_tensors =
      update_transformed_shape_values;

    if (update_flags & update_values)
      {
        values.resize(this->n_dofs_per_cell());
        data.shape_values.reinit(this->n_dofs_per_cell(), n_q_points);
        if (update_transformed_shape_values)
          data.transformed_shape_values.resize(n_q_points);
      }

    if (update_flags & update_gradients)
      {
        grads.resize(this->n_dofs_per_cell());
        data.shape_grads.reinit(this->n_dofs_per_cell(), n_q_points);
        data.transformed_shape_grads.resize(n_q_points);

        if (update_transformed_shape_grads)
          data.untransformed_shape_grads.resize(n_q_points);
      }

    if (update_flags & update_hessians)
      {
        grad_grads.resize(this->n_dofs_per_cell());
        data.shape_grad_grads.reinit(this->n_dofs_per_cell(), n_q_points);
        data.transformed_shape_hessians.resize(n_q_points);
        if (update_transformed_shape_hessian_tensors)
          data.untransformed_shape_hessian_tensors.resize(n_q_points);
      }

    // Compute shape function values
    // and derivatives and hessians on
    // the reference cell.
    // Make sure, that for the
    // node values N_i holds
    // N_i(v_j)=\delta_ij for all basis
    // functions v_j
    if (update_flags & (update_values | update_gradients))
      for (unsigned int k = 0; k < n_q_points; ++k)
        {
          poly_space->evaluate(quadrature.point(k),
                               values,
                               grads,
                               grad_grads,
                               third_derivatives,
                               fourth_derivatives);

          if (update_flags & update_values)
            {
              if (inverse_node_matrix.n_cols() == 0)
                for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
                  data.shape_values[i][k] = values[i];
              else
                for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
                  {
                    Tensor<1, dim> add_values;
                    for (unsigned int j = 0; j < this->n_dofs_per_cell(); ++j)
                      add_values += inverse_node_matrix(j, i) * values[j];
                    data.shape_values[i][k] = add_values;
                  }
            }

          if (update_flags & update_gradients)
            {
              if (inverse_node_matrix.n_cols() == 0)
                for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
                  data.shape_grads[i][k] = grads[i];
              else
                for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
                  {
                    Tensor<2, dim> add_grads;
                    for (unsigned int j = 0; j < this->n_dofs_per_cell(); ++j)
                      add_grads += inverse_node_matrix(j, i) * grads[j];
                    data.shape_grads[i][k] = add_grads;
                  }
            }

          if (update_flags & update_hessians)
            {
              if (inverse_node_matrix.n_cols() == 0)
                for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
                  data.shape_grad_grads[i][k] = grad_grads[i];
              else
                for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
                  {
                    Tensor<3, dim> add_grad_grads;
                    for (unsigned int j = 0; j < this->n_dofs_per_cell(); ++j)
                      add_grad_grads +=
                        inverse_node_matrix(j, i) * grad_grads[j];
                    data.shape_grad_grads[i][k] = add_grad_grads;
                  }
            }
        }
    return data_ptr;
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
   * Fields of cell-independent data for FE_PolyTensor. Stores the values of
   * the shape functions and their derivatives on the reference cell for later
   * use.
   *
   * All tables are organized in a way, that the value for shape function
   * <i>i</i> at quadrature point <i>k</i> is accessed by indices
   * <i>(i,k)</i>.
   */
  class InternalData : public FiniteElement<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * Array with shape function values in quadrature points. There is one row
     * for each shape function, containing values for each quadrature point.
     */
    Table<2, Tensor<1, dim>> shape_values;

    /**
     * Array with shape function gradients in quadrature points. There is one
     * row for each shape function, containing values for each quadrature
     * point.
     */
    Table<2, DerivativeForm<1, dim, spacedim>> shape_grads;

    /**
     * Array with shape function hessians in quadrature points. There is one
     * row for each shape function, containing values for each quadrature
     * point.
     */
    Table<2, DerivativeForm<2, dim, spacedim>> shape_grad_grads;

    /**
     * Scratch arrays for intermediate computations
     */
    mutable std::vector<double>              dof_sign_change;
    mutable std::vector<Tensor<1, spacedim>> transformed_shape_values;
    // for shape_gradient computations
    mutable std::vector<Tensor<2, spacedim>> transformed_shape_grads;
    mutable std::vector<Tensor<2, dim>>      untransformed_shape_grads;
    // for shape_hessian computations
    mutable std::vector<Tensor<3, spacedim>> transformed_shape_hessians;
    mutable std::vector<Tensor<3, dim>> untransformed_shape_hessian_tensors;
  };



  /**
   * A copy of the object passed to the constructor that describes the
   * polynomial space.
   */
  const std::unique_ptr<const TensorPolynomialsBase<dim>> poly_space;

  /**
   * The inverse of the matrix <i>a<sub>ij</sub></i> of node functionals
   * <i>N<sub>i</sub></i> applied to polynomial <i>p<sub>j</sub></i>. This
   * matrix is used to convert polynomials in the "raw" basis provided in
   * #poly_space to the basis dual to the node functionals on the reference
   * cell.
   *
   * This object is not filled by FE_PolyTensor, but is a chance for a derived
   * class to allow for reorganization of the basis functions. If it is left
   * empty, the basis in #poly_space is used, i.e., it is assumed that
   * #poly_space already satisfies the Kronecker delta property with
   * regard to the node functionals.
   */
  FullMatrix<double> inverse_node_matrix;

  /**
   * A mutex to be used to guard access to the variables below.
   */
  mutable Threads::Mutex cache_mutex;

  /**
   * If a shape function is computed at a single point, we must compute all of
   * them to apply #inverse_node_matrix. In order to avoid too much overhead,
   * we cache the point and the function values for the next evaluation.
   */
  mutable Point<dim> cached_point;

  /**
   * Cached shape function values after call to shape_value_component().
   */
  mutable std::vector<Tensor<1, dim>> cached_values;

  /**
   * Cached shape function gradients after call to shape_grad_component().
   */
  mutable std::vector<Tensor<2, dim>> cached_grads;

  /**
   * Cached second derivatives of shape functions after call to
   * shape_grad_grad_component().
   */
  mutable std::vector<Tensor<3, dim>> cached_grad_grads;
};

DEAL_II_NAMESPACE_CLOSE

#endif
