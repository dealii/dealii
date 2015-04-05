// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2015 by the deal.II authors
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

#ifndef __deal2__mapping_fe_h
#define __deal2__mapping_fe_h


#include <deal.II/base/config.h>
#include <deal.II/base/table.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/thread_management.h>


DEAL_II_NAMESPACE_OPEN


/*!@addtogroup mapping */
/*@{*/

/**
 * The MappingFE is a generalization of the MappingQEulerian class, for arbitrary
 * vectorial finite elements. The main difference is that this class uses a vector
 * of absolute positions instead of a vector of displacement.
 * In particular we think of a collections of a FE_Q or
 * Bezier finite element (FE_Bernstein) repeated a number of times equal to the space
 * dimension. The idea is to construct the mapping using a vector of control
 * points, a DoFHandler associated to the geometry of the problem and a
 * ComponentMask that tells us which components to use for the mapping.
 * This mapping will grab from the DoFHandler the finite element, or better
 * the collection of finite elements, to compute the mapping shape functions.
 * So we will have two different Finite Element and DoFHandler, one for the
 * solution field and one to describe the geometry of the problem. Historically
 * in the deal.II library there was not such a distinction. The differences
 * between this mapping and the MappingQ class are quite important.
 * The MappingFE, being a generalization, requires a higher level of abstraction.
 * This is the reason why it takes a DoFHandler and a vector of control points
 * that are the coefficients of the shape function (so in general it is a vector
 * of coefficient).
 *
 *
 * Typically, the DoFHandler operates on a finite element that is constructed
 * as a system element (FESystem) from continuous FE_Q() objects. An example
 * is shown below:
 * @code
 *    const FE_Q<dim,spacedim> feq(1);
 *    const FESystem<dim,spacedim> fesystem(feq, spacedim);
 *    DoFHandler<dim,spacedim> dhq(triangulation);
 *    dhq.distribute_dofs(fesystem);
 *    Vector<double> eulerq(dhq.n_dofs());
 *    const ComponentMask mask(spacedim, true);
 *    MappingFE<dim,spacedim> map(eulerq, dhq, mask);
 *    map.update_euler_vector_using_triangulation(eulerq);
 * @endcode


 *
 * @author Luca Heltai, Marco Tezzele 2013, 2015
 */
template <int dim, int spacedim=dim,
          class DH=DoFHandler<dim,spacedim>,
          class VECTOR=Vector<double> >
class MappingFE : public Mapping<dim,spacedim>
{
public:
  /**
   * Constructor. The first argument is a VECTOR that specifies the
   * transformation of the domain from the reference to the current
   * configuration. This is filled calling the method
   * update_euler_vector_using_triangulation.
   */
  MappingFE (const VECTOR  &euler_vector,
             const DH      &euler_dof_handler,
             const ComponentMask mask=ComponentMask());

  /**
   * Copy constructor. Performs a deep copy, i.e. duplicates what #tensor_pols
   * points to instead of simply copying the #tensor_pols pointer as done by a
   * default copy constructor.
   */
  MappingFE (const MappingFE<dim,spacedim,DH,VECTOR> &mapping);

  /**
   * Destructor.
   */
  virtual ~MappingFE ();

  /** Fill the euler vector with
  the information coming from
  the triangulation. Makes this
  map equivalent to MappingQ1,
  and it works ONLY if the
  underlying fe has support
  points. */
  void update_euler_vector_using_triangulation(VECTOR &vector);



  /**
   * Transforms the point @p p on the unit cell to the point @p p_real on the
   * real cell @p cell and returns @p p_real.
   */
  virtual Point<spacedim>
  transform_unit_to_real_cell (
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    const Point<dim>                                 &p) const;

  /**
   * Transforms the point @p p on the real cell to the point @p p_unit on the
   * unit cell @p cell and returns @p p_unit.
   *
   * Uses Newton iteration and the @p transform_unit_to_real_cell function.
   *
   * In the codimension one case, this function returns the normal projection
   * of the real point @p p on the curve or surface identified by the @p cell.
   */
  virtual Point<dim>
  transform_real_to_unit_cell (
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    const Point<spacedim>                            &p) const;


  virtual void
  transform (const VectorSlice<const std::vector<Tensor<1,dim> > > input,
             VectorSlice<std::vector<Tensor<1,spacedim> > > output,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const MappingType type) const;

  virtual void
  transform (const VectorSlice<const std::vector<DerivativeForm<1, dim, spacedim> > >    input,
             VectorSlice<std::vector<Tensor<2,spacedim> > > output,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const MappingType type) const;

  virtual
  void
  transform (const VectorSlice<const std::vector<Tensor<2, dim> > >     input,
             VectorSlice<std::vector<Tensor<2,spacedim> > >             output,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const MappingType type) const;



  /**
   * Return the degree of the mapping, i.e. the value which was passed to the
   * constructor.
   */
  unsigned int get_degree () const;

  /**
   * Return the ComponentMask of the mapping, i.e. which components to use for
   * the mapping.
   */
  ComponentMask get_fe_mask () const;

  /**
   * Return a pointer to a copy of the present object. The caller of this copy
   * then assumes ownership of it.
   */
  virtual
  Mapping<dim,spacedim> *clone () const;


  /**
   * Storage for internal data of
   * d-linear transformation.
   */
  class InternalData : public Mapping<dim,spacedim>::InternalDataBase
  {
  public:
    /**
     * Constructor.
     */
    InternalData(const FiniteElement<dim,spacedim> &fe,
                 const ComponentMask mask);

    /**
     * Shape function at quadrature
     * point. Shape functions are
     * in tensor product order, so
     * vertices must be reordered
     * to obtain transformation.
     */
    double shape (const unsigned int qpoint,
                  const unsigned int shape_nr) const;

    /**
     * Shape function at quadrature
     * point. See above.
     */
    double &shape (const unsigned int qpoint,
                   const unsigned int shape_nr);

    /**
     * Gradient of shape function
     * in quadrature point. See
     * above.
     */
    Tensor<1,dim> derivative (const unsigned int qpoint,
                              const unsigned int shape_nr) const;

    /**
     * Gradient of shape function
     * in quadrature point. See
     * above.
     */
    Tensor<1,dim> &derivative (const unsigned int qpoint,
                               const unsigned int shape_nr);

    /**
     * Second derivative of shape
     * function in quadrature
     * point. See above.
     */
    Tensor<2,dim> second_derivative (const unsigned int qpoint,
                                     const unsigned int shape_nr) const;

    /**
     * Second derivative of shape
     * function in quadrature
     * point. See above.
     */
    Tensor<2,dim> &second_derivative (const unsigned int qpoint,
                                      const unsigned int shape_nr);

    /**
     * Return an estimate (in
     * bytes) or the memory
     * consumption of this
     * object.
     */
    virtual std::size_t memory_consumption () const;

    /**
     * Values of shape
     * functions. Access by
     * function @p shape.
     *
     * Computed once.
     */
    std::vector<double> shape_values;

    /**
     * Values of shape function
     * derivatives. Access by
     * function @p derivative.
     *
     * Computed once.
     */
    std::vector<Tensor<1,dim> > shape_derivatives;

    /**
     * Values of shape function
     * second derivatives. Access
     * by function
     * @p second_derivative.
     *
     * Computed once.
     */
    std::vector<Tensor<2,dim> > shape_second_derivatives;

    /**
     * Tensors of covariant
     * transformation at each of
     * the quadrature points. The
     * matrix stored is the
     * Jacobian * G^{-1},
     * where G = Jacobian^{t} * Jacobian,
     * is the first fundamental
     * form of the map;
     * if dim=spacedim then
     * it reduces to the transpose of the
     * inverse of the Jacobian
     * matrix, which itself is
     * stored in the
     * @p contravariant field of
     * this structure.
     *
     * Computed on each cell.
     */
    std::vector<DerivativeForm<1,dim, spacedim > >  covariant;

    /**
     * Tensors of contravariant
     * transformation at each of
     * the quadrature points. The
     * contravariant matrix is
     * the Jacobian of the
     * transformation,
     * i.e. $J_{ij}=dx_i/d\hat x_j$.
     *
     * Computed on each cell.
     */
    std::vector< DerivativeForm<1,dim,spacedim> > contravariant;

    /**
     * Unit tangential vectors. Used
     * for the computation of
     * boundary forms and normal
     * vectors.
     *
     * This vector has
     * (dim-1)GeometryInfo::faces_per_cell
     * entries. The first
     * GeometryInfo::faces_per_cell
     * contain the vectors in the first
     * tangential direction for each
     * face; the second set of
     * GeometryInfo::faces_per_cell
     * entries contain the vectors in the
     * second tangential direction (only
     * in 3d, since there we have 2
     * tangential directions per face),
     * etc.
     *
     * Filled once.
     */
    std::vector<std::vector<Tensor<1,dim> > > unit_tangentials;

    /**
     * Auxiliary vectors for internal use.
     */
    std::vector<std::vector<Tensor<1,spacedim> > > aux;

    /**
     * Number of shape
     * functions. If this is a Q1
     * mapping, then it is simply
     * the number of vertices per
     * cell. However, since also
     * derived classes use this
     * class (e.g. the
     * Mapping_Q() class),
     * the number of shape
     * functions may also be
     * different.
     */
    unsigned int n_shape_functions;

    ComponentMask mask;
  };


  /**
   * Transforms a point @p p on
   * the unit cell to the point
   * @p p_real on the real cell
   * @p cell and returns @p p_real.
   *
   * This function is called by
   * @p transform_unit_to_real_cell
   * and multiple times (through the
   * Newton iteration) by
   * @p transform_real_to_unit_cell_internal.
   *
   * Takes a reference to an
   * @p InternalData that must
   * already include the shape
   * values at point @p p and the
   * mapping support points of the
   * cell.
   *
   * This @p InternalData argument
   * avoids multiple computations
   * of the shape values at point
   * @p p and especially multiple
   * computations of the mapping
   * support points.
   */
  Point<spacedim>
  transform_unit_to_real_cell_internal (const InternalData &mdata) const;


  /**
   * Transforms the point @p p on
   * the real cell to the corresponding
   * point on the unit cell
   * @p cell by a Newton
   * iteration.
   *
   * Takes a reference to an
   * @p InternalData that is
   * assumed to be previously
   * created by the @p get_data
   * function with @p UpdateFlags
   * including
   * @p update_transformation_values
   * and
   * @p update_transformation_gradients
   * and a one point Quadrature
   * that includes the given
   * initial guess for the
   * transformation
   * @p initial_p_unit.  Hence this
   * function assumes that
   * @p mdata already includes the
   * transformation shape values
   * and gradients computed at
   * @p initial_p_unit.
   *
   * @p mdata will be changed by
   * this function.
   */
  Point<dim>
  transform_real_to_unit_cell_internal (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                                        const Point<spacedim> &p,
                                        const Point<dim> &initial_p_unit,
                                        InternalData &mdata) const;

  /**
   * Do the computation for the
   * <tt>fill_*</tt> functions.
   */
  void compute_fill (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                     const unsigned int      npts,
                     const typename QProjector<dim>::DataSetDescriptor data_set,
                     const CellSimilarity::Similarity cell_similarity,
                     InternalData           &data,
                     std::vector<Point<spacedim> > &quadrature_points) const;


  /**
   * Do the computation for the
   * <tt>fill_*</tt> functions.
   */
  void compute_fill_face (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          const unsigned int      face_no,
                          const unsigned int      subface_no,
                          const unsigned int      npts,
                          const typename QProjector<dim>::DataSetDescriptor data_set,
                          const std::vector<double>   &weights,
                          InternalData           &mapping_data,
                          std::vector<Point<spacedim> >    &quadrature_points,
                          std::vector<double>         &JxW_values,
                          std::vector<Tensor<1,spacedim> > &boundary_form,
                          std::vector<Point<spacedim> > &normal_vectors,
                          std::vector<DerivativeForm<1,dim,spacedim> > &jacobians,
                          std::vector<DerivativeForm<1,spacedim,dim> > &inverse_jacobians) const;


  /**
   * Always returns @p false.
   */
  virtual
  bool preserves_vertex_locations () const;

  DeclException0(ExcInactiveCell);

protected:
  /**
   * Implementation of the interface in Mapping.
   */
  virtual void
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                  const Quadrature<dim>                                     &quadrature,
                  typename Mapping<dim,spacedim>::InternalDataBase          &mapping_data,
                  typename std::vector<Point<spacedim> >                    &quadrature_points,
                  std::vector<double>                                       &JxW_values,
                  std::vector<DerivativeForm<1,dim,spacedim> >       &jacobians,
                  std::vector<DerivativeForm<2,dim,spacedim> >       &jacobian_grads,
                  std::vector<DerivativeForm<1,spacedim,dim> >      &inverse_jacobians,
                  std::vector<Point<spacedim> >                             &cell_normal_vectors,
                  CellSimilarity::Similarity                           &cell_similarity) const ;

  /**
   * Implementation of the interface in Mapping.
   */
  virtual void
  fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                       const unsigned int face_no,
                       const Quadrature<dim-1>& quadrature,
                       typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
                       typename std::vector<Point<spacedim> >        &quadrature_points,
                       std::vector<double>             &JxW_values,
                       std::vector<Tensor<1,spacedim> >             &exterior_forms,
                       std::vector<Point<spacedim> >                &normal_vectors,
                       std::vector<DerivativeForm<1,dim,spacedim> > &jacobians,
                       std::vector<DerivativeForm<1,spacedim,dim> > &inverse_jacobians) const ;

  /**
   * Implementation of the interface in Mapping.
   */
  virtual void
  fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          const unsigned int face_no,
                          const unsigned int sub_no,
                          const Quadrature<dim-1>& quadrature,
                          typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
                          typename std::vector<Point<spacedim> >        &quadrature_points,
                          std::vector<double>             &JxW_values,
                          std::vector<Tensor<1,spacedim> > &exterior_forms,
                          std::vector<Point<spacedim> >    &normal_vectors,
                          std::vector<DerivativeForm<1,dim,spacedim> > &jacobians,
                          std::vector<DerivativeForm<1,spacedim,dim> > &inverse_jacobians) const ;


  /**
     This function and the next allow to generate the transform require by
     the virtual transform() in mapping, but unfortunately in C++ one cannot
     declare a virtual template function.
  */
  template < int rank >
  void
  transform_fields(const VectorSlice<const std::vector<Tensor<rank,dim>      > > input,
                   VectorSlice<      std::vector<Tensor<rank,spacedim> > > output,
                   const typename Mapping<dim,spacedim>::InternalDataBase &internal,
                   const MappingType type) const;


  /**
       see doc in transform_fields
   */
  template < int rank >
  void
  transform_differential_forms(
    const VectorSlice<const std::vector<DerivativeForm<rank, dim,spacedim> > >    input,
    VectorSlice<std::vector<DerivativeForm<rank, spacedim,spacedim> > > output,
    const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
    const MappingType mapping_type) const;


protected:
  /**
  * Reference to the vector of shifts.
  */

  SmartPointer<const VECTOR, MappingFE<dim,spacedim,DH,VECTOR> >euler_vector;
  /**
   * A FiniteElement object which is only needed in 3D, since it knows how to reorder
   * shape functions/DoFs on non-standard faces. This is used to reorder
   * support points in the same way. We could make this a pointer to prevent
   * construction in 1D and 2D, but since memory and time requirements are not
   * particularly high this seems unnecessary at the moment.
   */
  SmartPointer<const FiniteElement<dim,spacedim>, MappingFE<dim,spacedim,DH,VECTOR> > fe;


  /**
   * Pointer to the DoFHandler to which the mapping vector is associated.
   */
  SmartPointer<const DH,MappingFE<dim,spacedim,DH,VECTOR> >euler_dof_handler;



private:
//
  /**
   * Update internal degrees of
   * freedom. */
  void update_internal_dofs(const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;


  mutable std::vector<double> local_dofs;

  mutable std::vector<unsigned int> dof_indices;

  /**
   * Mutex to protect local_dofs.
   */

  mutable Threads::Mutex mutex;


  virtual void
  compute_shapes_virtual (const std::vector<Point<dim> > &unit_points,
                          typename MappingFE<dim, spacedim>::InternalData &data) const;

  UpdateFlags
  update_once (const UpdateFlags in) const;

  UpdateFlags
  update_each (const UpdateFlags in) const;

  void
  compute_data (const UpdateFlags      update_flags,
                const Quadrature<dim>  &q,
                const unsigned int     n_original_q_points,
                InternalData           &data) const;

  void
  compute_face_data (const UpdateFlags      update_flags,
                     const Quadrature<dim>  &q,
                     const unsigned int     n_original_q_points,
                     InternalData           &data) const;

  virtual
  typename Mapping<dim,spacedim>::InternalDataBase *
  get_data (const UpdateFlags,
            const Quadrature<dim> &quadrature) const;

  virtual
  typename Mapping<dim,spacedim>::InternalDataBase *
  get_face_data (const UpdateFlags flags,
                 const Quadrature<dim-1>& quadrature) const;

  virtual
  typename Mapping<dim,spacedim>::InternalDataBase *
  get_subface_data (const UpdateFlags flags,
                    const Quadrature<dim-1>& quadrature) const;


  /*
   * Which components to use for the mapping.
   */
  const ComponentMask fe_mask;


  /**
   * Mapping between indices in the FE space and the real space. This vector contains one
   * index for each component of the finite element space. If the index is one for which
   * the ComponentMask which is used to construct this element is false, then
   * numbers::invalid_unsigned_int is returned, otherwise the component in real space is
   * returned. For example, if we construct the mapping using ComponentMask(spacedim, true),
   * then this vector contains {0,1,2} in spacedim = 3.
   */
  std::vector<unsigned int> fe_to_real;


  /**
   * Declare other MappingFE classes friends.
   */
  template <int,int,class,class> friend class MappingFE;
};

/*@}*/

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN


template<int dim, int spacedim, class DH, class VECTOR>
inline
double
MappingFE<dim,spacedim,DH,VECTOR>::InternalData::shape (const unsigned int qpoint,
                                                        const unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_values.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_values.size()));
  return shape_values [qpoint*n_shape_functions + shape_nr];
}



template<int dim, int spacedim, class DH, class VECTOR>
inline
double &
MappingFE<dim,spacedim,DH,VECTOR>::InternalData::shape (const unsigned int qpoint,
                                                        const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_values.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_values.size()));
  return shape_values [qpoint*n_shape_functions + shape_nr];
}


template<int dim, int spacedim, class DH, class VECTOR>
inline
Tensor<1,dim>
MappingFE<dim,spacedim,DH,VECTOR>::InternalData::derivative (const unsigned int qpoint,
    const unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_derivatives.size()));
  return shape_derivatives [qpoint*n_shape_functions + shape_nr];
}



template<int dim, int spacedim, class DH, class VECTOR>
inline
Tensor<1,dim> &
MappingFE<dim,spacedim,DH,VECTOR>::InternalData::derivative (const unsigned int qpoint,
    const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_derivatives.size()));
  return shape_derivatives [qpoint*n_shape_functions + shape_nr];
}


template <int dim, int spacedim, class DH, class VECTOR>
inline
Tensor<2,dim>
MappingFE<dim,spacedim,DH,VECTOR>::InternalData::second_derivative (const unsigned int qpoint,
    const unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_second_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_second_derivatives.size()));
  return shape_second_derivatives [qpoint*n_shape_functions + shape_nr];
}



template <int dim, int spacedim, class DH, class VECTOR>
inline
Tensor<2,dim> &
MappingFE<dim,spacedim,DH,VECTOR>::InternalData::second_derivative (const unsigned int qpoint,
    const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_second_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_second_derivatives.size()));
  return shape_second_derivatives [qpoint*n_shape_functions + shape_nr];
}


template <int dim, int spacedim, class DH, class VECTOR>
inline
bool
MappingFE<dim,spacedim,DH,VECTOR>::preserves_vertex_locations () const
{
  return false;
}




#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
