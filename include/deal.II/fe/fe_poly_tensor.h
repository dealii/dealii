// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2016 by the deal.II authors
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

#ifndef dealii__fe_poly_tensor_h
#define dealii__fe_poly_tensor_h


#include <deal.II/base/config.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/fe/fe.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/thread_management.h>


DEAL_II_NAMESPACE_OPEN

/**
 * This class gives a unified framework for the implementation of
 * FiniteElement classes based on Tensor valued polynomial spaces like
 * PolynomialsBDM and PolynomialsRaviartThomas.
 *
 * Every class that implements following function can be used as template
 * parameter PolynomialType.
 *
 * @code
 * void compute (const Point<dim>            &unit_point,
 *               std::vector<Tensor<1,dim> > &values,
 *               std::vector<Tensor<2,dim> > &grads,
 *               std::vector<Tensor<3,dim> > &grad_grads) const;
 * @endcode
 *
 * In many cases, the node functionals depend on the shape of the mesh cell,
 * since they evaluate normal or tangential components on the faces. In order
 * to allow for a set of transformations, the variable #mapping_type has been
 * introduced. It should also be set in the constructor of a derived class.
 *
 * This class is not a fully implemented FiniteElement class, but implements
 * some common features of vector valued elements based on vector valued
 * polynomial classes. What's missing here in particular is information on the
 * topological location of the node values.
 *
 * For more information on the template parameter <tt>spacedim</tt>, see the
 * documentation for the class Triangulation.
 *
 * <h3>Deriving classes</h3>
 *
 * Any derived class must decide on the polynomial space to use.  This
 * polynomial space should be implemented simply as a set of vector valued
 * polynomials like PolynomialsBDM and PolynomialsRaviartThomas.  In order to
 * facilitate this implementation, the basis of this space may be arbitrary.
 *
 * <h4>Determining the correct basis</h4>
 *
 * In most cases, the set of desired node values $N_i$ and the basis functions
 * $v_j$ will not fulfill the interpolation condition $N_i(v_j) =
 * \delta_{ij}$.
 *
 * The use of the member data #inverse_node_matrix allows to compute the basis
 * $v_j$ automatically, after the node values for each original basis function
 * of the polynomial space have been computed.
 *
 * Therefore, the constructor of a derived class should have a structure like
 * this (example for interpolation in support points):
 *
 * @code
 *  fill_support_points();
 *
 *  const unsigned int n_dofs = this->dofs_per_cell;
 *  FullMatrix<double> N(n_dofs, n_dofs);
 *
 *  for (unsigned int i=0;i<n_dofs;++i)
 *    {
 *      const Point<dim>& p = this->unit_support_point[i];
 *
 *      for (unsigned int j=0;j<n_dofs;++j)
 *        for (unsigned int d=0;d<dim;++d)
 *          N(i,j) += node_vector[i][d]
 *                  * this->shape_value_component(j, p, d);
 *    }
 *
 *  this->inverse_node_matrix.reinit(n_dofs, n_dofs);
 *  this->inverse_node_matrix.invert(N);
 * @endcode
 *
 * @note The matrix #inverse_node_matrix should have dimensions zero before
 * this piece of code is executed. Only then, shape_value_component() will
 * return the raw polynomial <i>j</i> as defined in the polynomial space
 * PolynomialType.
 *
 * <h4>Setting the transformation</h4>
 *
 * In most cases, vector valued basis functions must be transformed when
 * mapped from the reference cell to the actual grid cell. These
 * transformations can be selected from the set MappingType and stored in
 * #mapping_type. Therefore, each constructor should contain a line like:
 * @code
 * this->mapping_type = this->mapping_none;
 * @endcode
 *
 * @see PolynomialsBDM, PolynomialsRaviartThomas
 * @ingroup febase
 * @author Guido Kanschat
 * @date 2005
 */
template <class PolynomialType, int dim, int spacedim=dim>
class FE_PolyTensor : public FiniteElement<dim,spacedim>
{
public:
  /**
   * Constructor.
   *
   * @arg @c degree: constructor argument for poly. May be different from @p
   * fe_data.degree.
   */
  FE_PolyTensor (const unsigned int degree,
                 const FiniteElementData<dim> &fe_data,
                 const std::vector<bool> &restriction_is_additive_flags,
                 const std::vector<ComponentMask> &nonzero_components);

  // for documentation, see the FiniteElement base class
  virtual
  UpdateFlags
  requires_update_flags (const UpdateFlags update_flags) const;

  /**
   * Compute the (scalar) value of shape function @p i at the given quadrature
   * point @p p. Since the elements represented by this class are vector
   * valued, there is no such scalar value and the function therefore throws
   * an exception.
   */
  virtual double shape_value (const unsigned int i,
                              const Point<dim> &p) const;

  // documentation inherited from the base class
  virtual double shape_value_component (const unsigned int i,
                                        const Point<dim> &p,
                                        const unsigned int component) const;

  /**
   * Compute the gradient of (scalar) shape function @p i at the given
   * quadrature point @p p. Since the elements represented by this class are
   * vector valued, there is no such scalar value and the function therefore
   * throws an exception.
   */
  virtual Tensor<1,dim> shape_grad (const unsigned int  i,
                                    const Point<dim>   &p) const;

  // documentation inherited from the base class
  virtual Tensor<1,dim> shape_grad_component (const unsigned int i,
                                              const Point<dim> &p,
                                              const unsigned int component) const;

  /**
   * Compute the Hessian of (scalar) shape function @p i at the given
   * quadrature point @p p. Since the elements represented by this class are
   * vector valued, there is no such scalar value and the function therefore
   * throws an exception.
   */
  virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
                                         const Point<dim> &p) const;

  // documentation inherited from the base class
  virtual Tensor<2,dim> shape_grad_grad_component (const unsigned int i,
                                                   const Point<dim> &p,
                                                   const unsigned int component) const;

protected:
  /**
   * The mapping type to be used to map shape functions from the reference
   * cell to the mesh cell.
   */
  MappingType mapping_type;


  /* NOTE: The following function has its definition inlined into the class declaration
     because we otherwise run into a compiler error with MS Visual Studio. */
  virtual
  typename FiniteElement<dim,spacedim>::InternalDataBase *
  get_data(const UpdateFlags                                                    update_flags,
           const Mapping<dim,spacedim>                                         &/*mapping*/,
           const Quadrature<dim>                                               &quadrature,
           dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &/*output_data*/) const
  {
    // generate a new data object and
    // initialize some fields
    InternalData *data = new InternalData;
    data->update_each = requires_update_flags(update_flags);

    const unsigned int n_q_points = quadrature.size();

    // some scratch arrays
    std::vector<Tensor<1,dim> > values(0);
    std::vector<Tensor<2,dim> > grads(0);
    std::vector<Tensor<3,dim> > grad_grads(0);
    std::vector<Tensor<4,dim> > third_derivatives(0);
    std::vector<Tensor<5,dim> > fourth_derivatives(0);

    if (update_flags & (update_values | update_gradients | update_hessians) )
      data->sign_change.resize (this->dofs_per_cell);

    // initialize fields only if really
    // necessary. otherwise, don't
    // allocate memory
    if (update_flags & update_values)
      {
        values.resize (this->dofs_per_cell);
        data->shape_values.reinit (this->dofs_per_cell, n_q_points);
        if (mapping_type != mapping_none)
          data->transformed_shape_values.resize (n_q_points);
      }

    if (update_flags & update_gradients)
      {
        grads.resize (this->dofs_per_cell);
        data->shape_grads.reinit (this->dofs_per_cell, n_q_points);
        data->transformed_shape_grads.resize (n_q_points);

        if ( (mapping_type == mapping_raviart_thomas)
             ||
             (mapping_type == mapping_piola)
             ||
             (mapping_type == mapping_nedelec)
             ||
             (mapping_type == mapping_contravariant))
          data->untransformed_shape_grads.resize(n_q_points);
      }

    if (update_flags & update_hessians)
      {
        grad_grads.resize (this->dofs_per_cell);
        data->shape_grad_grads.reinit (this->dofs_per_cell, n_q_points);
        data->transformed_shape_hessians.resize (n_q_points);
        if ( mapping_type != mapping_none )
          data->untransformed_shape_hessian_tensors.resize(n_q_points);
      }

    // Compute shape function values
    // and derivatives and hessians on
    // the reference cell.
    // Make sure, that for the
    // node values N_i holds
    // N_i(v_j)=\delta_ij for all basis
    // functions v_j
    if (update_flags & (update_values | update_gradients))
      for (unsigned int k=0; k<n_q_points; ++k)
        {
          poly_space.compute(quadrature.point(k),
                             values, grads, grad_grads,
                             third_derivatives,
                             fourth_derivatives);

          if (update_flags & update_values)
            {
              if (inverse_node_matrix.n_cols() == 0)
                for (unsigned int i=0; i<this->dofs_per_cell; ++i)
                  data->shape_values[i][k] = values[i];
              else
                for (unsigned int i=0; i<this->dofs_per_cell; ++i)
                  {
                    Tensor<1,dim> add_values;
                    for (unsigned int j=0; j<this->dofs_per_cell; ++j)
                      add_values += inverse_node_matrix(j,i) * values[j];
                    data->shape_values[i][k] = add_values;
                  }
            }

          if (update_flags & update_gradients)
            {
              if (inverse_node_matrix.n_cols() == 0)
                for (unsigned int i=0; i<this->dofs_per_cell; ++i)
                  data->shape_grads[i][k] = grads[i];
              else
                for (unsigned int i=0; i<this->dofs_per_cell; ++i)
                  {
                    Tensor<2,dim> add_grads;
                    for (unsigned int j=0; j<this->dofs_per_cell; ++j)
                      add_grads += inverse_node_matrix(j,i) * grads[j];
                    data->shape_grads[i][k] = add_grads;
                  }
            }

          if (update_flags & update_hessians)
            {
              if (inverse_node_matrix.n_cols() == 0)
                for (unsigned int i=0; i<this->dofs_per_cell; ++i)
                  data->shape_grad_grads[i][k] = grad_grads[i];
              else
                for (unsigned int i=0; i<this->dofs_per_cell; ++i)
                  {
                    Tensor<3,dim> add_grad_grads;
                    for (unsigned int j=0; j<this->dofs_per_cell; ++j)
                      add_grad_grads += inverse_node_matrix(j,i) * grad_grads[j];
                    data->shape_grad_grads[i][k] = add_grad_grads;
                  }
            }

        }
    return data;
  }

  virtual
  void
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator           &cell,
                  const CellSimilarity::Similarity                                     cell_similarity,
                  const Quadrature<dim>                                               &quadrature,
                  const Mapping<dim,spacedim>                                         &mapping,
                  const typename Mapping<dim,spacedim>::InternalDataBase              &mapping_internal,
                  const dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &mapping_data,
                  const typename FiniteElement<dim,spacedim>::InternalDataBase        &fe_internal,
                  dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &output_data) const;

  virtual
  void
  fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator           &cell,
                       const unsigned int                                                   face_no,
                       const Quadrature<dim-1>                                             &quadrature,
                       const Mapping<dim,spacedim>                                         &mapping,
                       const typename Mapping<dim,spacedim>::InternalDataBase              &mapping_internal,
                       const dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &mapping_data,
                       const typename FiniteElement<dim,spacedim>::InternalDataBase        &fe_internal,
                       dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &output_data) const;

  virtual
  void
  fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator           &cell,
                          const unsigned int                                                   face_no,
                          const unsigned int                                                   sub_no,
                          const Quadrature<dim-1>                                             &quadrature,
                          const Mapping<dim,spacedim>                                         &mapping,
                          const typename Mapping<dim,spacedim>::InternalDataBase              &mapping_internal,
                          const dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &mapping_data,
                          const typename FiniteElement<dim,spacedim>::InternalDataBase        &fe_internal,
                          dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &output_data) const;

  /**
   * Fields of cell-independent data for FE_PolyTensor. Stores the values of
   * the shape functions and their derivatives on the reference cell for later
   * use.
   *
   * All tables are organized in a way, that the value for shape function
   * <i>i</i> at quadrature point <i>k</i> is accessed by indices
   * <i>(i,k)</i>.
   */
  class InternalData : public FiniteElement<dim,spacedim>::InternalDataBase
  {
  public:
    /**
     * Array with shape function values in quadrature points. There is one row
     * for each shape function, containing values for each quadrature point.
     */
    Table<2,Tensor<1,dim> > shape_values;

    /**
     * Array with shape function gradients in quadrature points. There is one
     * row for each shape function, containing values for each quadrature
     * point.
     */
    Table<2,DerivativeForm<1, dim, spacedim> > shape_grads;

    /**
     * Array with shape function hessians in quadrature points. There is one
     * row for each shape function, containing values for each quadrature
     * point.
     */
    Table<2,DerivativeForm<2, dim, spacedim> > shape_grad_grads;

    /**
     * Scratch arrays for intermediate computations
     */
    mutable std::vector<double>                sign_change;
    mutable std::vector<Tensor<1, spacedim> >  transformed_shape_values;
    // for shape_gradient computations
    mutable std::vector<Tensor<2, spacedim > > transformed_shape_grads;
    mutable std::vector<Tensor<2, dim > >      untransformed_shape_grads;
    // for shape_hessian computations
    mutable std::vector<Tensor<3, spacedim > > transformed_shape_hessians;
    mutable std::vector<Tensor<3, dim > >      untransformed_shape_hessian_tensors;
  };



  /**
   * The polynomial space. Its type is given by the template parameter
   * PolynomialType.
   */
  PolynomialType poly_space;

  /**
   * The inverse of the matrix <i>a<sub>ij</sub></i> of node values
   * <i>N<sub>i</sub></i> applied to polynomial <i>p<sub>j</sub></i>. This
   * matrix is used to convert polynomials in the "raw" basis provided in
   * #poly_space to the basis dual to the node functionals on the reference
   * cell.
   *
   * This object is not filled by FE_PolyTensor, but is a chance for a derived
   * class to allow for reorganization of the basis functions. If it is left
   * empty, the basis in #poly_space is used.
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
  mutable std::vector<Tensor<1,dim> > cached_values;

  /**
   * Cached shape function gradients after call to shape_grad_component().
   */
  mutable std::vector<Tensor<2,dim> > cached_grads;

  /**
   * Cached second derivatives of shape functions after call to
   * shape_grad_grad_component().
   */
  mutable std::vector<Tensor<3,dim> > cached_grad_grads;
};

DEAL_II_NAMESPACE_CLOSE

#endif
