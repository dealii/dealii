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

#ifndef dealii__mapping_fe_h
#define dealii__mapping_fe_h


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
 * The MappingFEField is a generalization of the MappingQEulerian class, for
 * arbitrary vector finite elements. The two main differences are that this
 * class uses a vector of absolute positions instead of a vector of
 * displacements, and it allows for arbitrary FiniteElement types, instead of
 * only FE_Q.
 *
 * This class effectively decouples the topology from the geometry, by
 * relegating all geometrical information to some components of a
 * FiniteElement vector field. The components that are used for the geometry
 * can be arbitrarily selected at construction time.
 *
 * The idea is to consider the Triangulation as a parameter configuration
 * space, on which we  construct an arbitrary geometrical mapping, using the
 * instruments of the deal.II library: a vector of degrees of freedom, a
 * DoFHandler associated to the geometry of the problem and a ComponentMask
 * that tells us which components of the FiniteElement to use for the mapping.
 *
 * Typically, the DoFHandler operates on a finite element that is constructed
 * as a system element (FESystem()) from continuous FE_Q() (for iso-parametric
 * discretizations) or FE_Bernstein() (for iso-geometric discretizations)
 * objects. An example is shown below:
 *
 * @code
 *    const FE_Q<dim,spacedim> feq(1);
 *    const FESystem<dim,spacedim> fesystem(feq, spacedim);
 *    DoFHandler<dim,spacedim> dhq(triangulation);
 *    dhq.distribute_dofs(fesystem);
 *    const ComponentMask mask(spacedim, true);
 *    Vector<double> eulerq(dhq.n_dofs());
 *    // Fills the euler vector with information from the Triangulation
 *    VectorTools::get_position_vector(dhq, eulerq, mask);
 *    MappingFEField<dim,spacedim> map(dhq, eulerq, mask);
 * @endcode
 *
 * @author Luca Heltai, Marco Tezzele 2013, 2015
 */
template <int dim, int spacedim=dim,
          class VECTOR=Vector<double>,
          class DH=DoFHandler<dim,spacedim> >
class MappingFEField : public Mapping<dim,spacedim>
{
public:
  /**
   * Constructor. The first argument is a VECTOR that specifies the
   * transformation of the domain from the reference to the current
   * configuration.
   *
   * In general this class decouples geometry from topology, allowing users to
   * define geometries which are only topologically equivalent to the
   * underlying Triangulation, but which may otherwise be arbitrary.
   * Differently from what happens in MappingQEulerian, the FiniteElement
   * field which is passed to the constructor is interpreted as an absolute
   * geometrical configuration, therefore one has to make sure that the
   * euler_vector actually represents a valid geometry (i.e., one with no
   * inverted cells, or with no zero-volume cells).
   *
   * If the underlying FiniteElement is a system of FE_Q(), and euler_vector
   * is initialized using VectorTools::get_position_vector(), then this class
   * is in all respects identical to MappingQ().
   *
   * The optional ComponentMask argument can be used to specify what
   * components of the FiniteElement to use for the geometrical
   * transformation. If no mask is specified at construction time, then a
   * default one is used, which makes this class works in the same way of
   * MappingQEulerian(), i.e., the first spacedim components of the
   * FiniteElement are assumed to represent the geometry of the problem.
   *
   * Notice that if a mask is specified, it has to match in size the
   * underlying FiniteElement, and it has to have exactly spacedim non-zero
   * elements, indicating the components (in order) of the FiniteElement which
   * will be used for the geometry.
   *
   * If an incompatible mask is passed, an exception is thrown.
   */
  MappingFEField (const DH      &euler_dof_handler,
                  const VECTOR  &euler_vector,
                  const ComponentMask mask=ComponentMask());

  /**
   * Copy constructor.
   */
  MappingFEField (const MappingFEField<dim,spacedim,VECTOR,DH> &mapping);


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

  /**
   * Reimplemented from Mapping. See the documentation of the base class for
   * detailed information.
   */
  virtual void
  transform (const VectorSlice<const std::vector<Tensor<1,dim> > > input,
             VectorSlice<std::vector<Tensor<1,spacedim> > > output,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const MappingType type) const;

  /**
   * Reimplemented from Mapping. See the documentation of the base class for
   * detailed information.
   */
  virtual void
  transform (const VectorSlice<const std::vector<DerivativeForm<1, dim, spacedim> > >    input,
             VectorSlice<std::vector<Tensor<2,spacedim> > > output,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const MappingType type) const;

  /**
   * Reimplemented from Mapping. See the documentation of the base class for
   * detailed information.
   */
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
  ComponentMask get_component_mask () const;

  /**
   * Return a pointer to a copy of the present object. The caller of this copy
   * then assumes ownership of it.
   */
  virtual
  Mapping<dim,spacedim> *clone () const;


  /**
   * Storage for internal data of this mapping. See Mapping::InternalDataBase
   * for an extensive description.
   *
   * This includes data that is computed once when the object is created
   * (in get_data()) as well as data the class wants to store from between
   * the call to fill_fe_values(), fill_fe_face_values(), or
   * fill_fe_subface_values() until possible later calls from the finite
   * element to functions such as transform(). The latter class of
   * member variables are marked as 'mutable', along with scratch arrays.
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
     * Shape function at quadrature point. Shape functions are in tensor
     * product order, so vertices must be reordered to obtain transformation.
     */
    const double &shape (const unsigned int qpoint,
                         const unsigned int shape_nr) const;

    /**
     * Shape function at quadrature point. See above.
     */
    double &shape (const unsigned int qpoint,
                   const unsigned int shape_nr);

    /**
     * Gradient of shape function in quadrature point. See above.
     */
    const Tensor<1,dim> &derivative (const unsigned int qpoint,
                                     const unsigned int shape_nr) const;

    /**
     * Gradient of shape function in quadrature point. See above.
     */
    Tensor<1,dim> &derivative (const unsigned int qpoint,
                               const unsigned int shape_nr);

    /**
     * Second derivative of shape function in quadrature point. See above.
     */
    const Tensor<2,dim> &second_derivative (const unsigned int qpoint,
                                            const unsigned int shape_nr) const;

    /**
     * Second derivative of shape function in quadrature point. See above.
     */
    Tensor<2,dim> &second_derivative (const unsigned int qpoint,
                                      const unsigned int shape_nr);

    /**
     * Return an estimate (in bytes) or the memory consumption of this object.
     */
    virtual std::size_t memory_consumption () const;

    /**
     * Values of shape functions. Access by function @p shape.
     *
     * Computed once.
     */
    std::vector<double> shape_values;

    /**
     * Values of shape function derivatives. Access by function @p derivative.
     *
     * Computed once.
     */
    std::vector<Tensor<1,dim> > shape_derivatives;

    /**
     * Values of shape function second derivatives. Access by function @p
     * second_derivative.
     *
     * Computed once.
     */
    std::vector<Tensor<2,dim> > shape_second_derivatives;

    /**
     * Unit tangential vectors. Used for the computation of boundary forms and
     * normal vectors.
     *
     * This vector has (dim-1)GeometryInfo::faces_per_cell entries. The first
     * GeometryInfo::faces_per_cell contain the vectors in the first
     * tangential direction for each face; the second set of
     * GeometryInfo::faces_per_cell entries contain the vectors in the second
     * tangential direction (only in 3d, since there we have 2 tangential
     * directions per face), etc.
     *
     * Filled once.
     */
    std::vector<std::vector<Tensor<1,dim> > > unit_tangentials;

    /**
     * Number of shape functions. If this is a Q1 mapping, then it is simply
     * the number of vertices per cell. However, since also derived classes
     * use this class (e.g. the Mapping_Q() class), the number of shape
     * functions may also be different.
     */
    unsigned int n_shape_functions;

    /**
     * Stores the mask given at construction time. If no mask was specified at
     * construction time, then a default one is used, which makes this class
     * works in the same way of MappingQEulerian(), i.e., the first spacedim
     * components of the FiniteElement are used for the euler_vector and the
     * euler_dh.
     *
     * If a mask is specified, then it has to match the underlying
     * FiniteElement, and it has to have exactly spacedim non-zero elements,
     * indicating the components (in order) of the FiniteElement which will be
     * used for the euler vector and the euler dof handler.
     */
    ComponentMask mask;

    /**
     * Tensors of covariant transformation at each of the quadrature points.
     * The matrix stored is the Jacobian * G^{-1}, where G = Jacobian^{t} *
     * Jacobian, is the first fundamental form of the map; if dim=spacedim
     * then it reduces to the transpose of the inverse of the Jacobian matrix,
     * which itself is stored in the @p contravariant field of this structure.
     *
     * Computed on each cell.
     */
    mutable std::vector<DerivativeForm<1,dim, spacedim > >  covariant;

    /**
     * Tensors of contravariant transformation at each of the quadrature
     * points. The contravariant matrix is the Jacobian of the transformation,
     * i.e. $J_{ij}=dx_i/d\hat x_j$.
     *
     * Computed on each cell.
     */
    mutable std::vector< DerivativeForm<1,dim,spacedim> > contravariant;

    /**
     * The determinant of the Jacobian in each quadrature point. Filled if
     * #update_volume_elements.
     */
    mutable std::vector<double> volume_elements;

    /**
     * Auxiliary vectors for internal use.
     */
    mutable std::vector<std::vector<Tensor<1,spacedim> > > aux;

    /**
     * Storage for the indices of the local degrees of freedom.
     */
    mutable std::vector<types::global_dof_index> local_dof_indices;

    /**
     * Storage for local degrees of freedom.
     */
    mutable std::vector<double> local_dof_values;
  };


  /**
   * Transforms a point @p p on the unit cell to the point @p p_real on the
   * real cell @p cell and returns @p p_real.
   *
   * This function is called by @p transform_unit_to_real_cell and multiple
   * times (through the Newton iteration) by @p
   * transform_real_to_unit_cell_internal.
   *
   * Takes a reference to an @p InternalData that must already include the
   * shape values at point @p p and the mapping support points of the cell.
   *
   * This @p InternalData argument avoids multiple computations of the shape
   * values at point @p p and especially multiple computations of the mapping
   * support points.
   */
  Point<spacedim>
  transform_unit_to_real_cell_internal (const InternalData &mdata) const;


  /**
   * Transforms the point @p p on the real cell to the corresponding point on
   * the unit cell @p cell by a Newton iteration.
   *
   * Takes a reference to an @p InternalData that is assumed to be previously
   * created by the @p get_data function with @p UpdateFlags including @p
   * update_transformation_values and @p update_transformation_gradients and a
   * one point Quadrature that includes the given initial guess for the
   * transformation @p initial_p_unit.  Hence this function assumes that @p
   * mdata already includes the transformation shape values and gradients
   * computed at @p initial_p_unit.
   *
   * @p mdata will be changed by this function.
   */
  Point<dim>
  transform_real_to_unit_cell_internal (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                                        const Point<spacedim> &p,
                                        const Point<dim> &initial_p_unit,
                                        InternalData &mdata) const;

  /**
   * Do the computation for the <tt>fill_*</tt> functions.
   */
  void compute_fill (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                     const unsigned int      npts,
                     const typename QProjector<dim>::DataSetDescriptor data_set,
                     const CellSimilarity::Similarity cell_similarity,
                     const InternalData           &data,
                     std::vector<Point<spacedim> > &quadrature_points) const;


  /**
   * Do the computation for the <tt>fill_*</tt> functions.
   */
  void compute_fill_face (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          const unsigned int      face_no,
                          const unsigned int      subface_no,
                          const unsigned int      npts,
                          const typename QProjector<dim>::DataSetDescriptor data_set,
                          const std::vector<double>   &weights,
                          const InternalData           &mapping_data,
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
   * Compute mapping-related information for a cell.
   * See the documentation of Mapping::fill_fe_values() for
   * a discussion of purpose, arguments, and return value of this function.
   */
  virtual
  CellSimilarity::Similarity
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                  const CellSimilarity::Similarity                           cell_similarity,
                  const Quadrature<dim>                                     &quadrature,
                  const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                  FEValuesData<dim,spacedim>                                &output_data) const;

  /**
   * Compute mapping-related information for a face of a cell.
   * See the documentation of Mapping::fill_fe_face_values() for
   * a discussion of purpose and arguments of this function.
   */
  virtual void
  fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                       const unsigned int                                         face_no,
                       const Quadrature<dim-1>                                   &quadrature,
                       const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                       FEValuesData<dim,spacedim>                                &output_data) const;

  /**
   * Compute mapping-related information for a child of a face of a cell.
   * See the documentation of Mapping::fill_fe_subface_values() for
   * a discussion of purpose and arguments of this function.
   */
  virtual void
  fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          const unsigned int                                         face_no,
                          const unsigned int                                         subface_no,
                          const Quadrature<dim-1>                                   &quadrature,
                          const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                          FEValuesData<dim,spacedim>                                &output_data) const;

  /**
   * This function and the next allow to generate the transform require by the
   * virtual transform() in mapping, but unfortunately in C++ one cannot
   * declare a virtual template function.
   */
  template < int rank >
  void
  transform_fields(const VectorSlice<const std::vector<Tensor<rank,dim>      > > input,
                   VectorSlice<      std::vector<Tensor<rank,spacedim> > > output,
                   const typename Mapping<dim,spacedim>::InternalDataBase &internal,
                   const MappingType type) const;


  /**
   * see doc in transform_fields
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

  SmartPointer<const VECTOR, MappingFEField<dim,spacedim,DH,VECTOR> >euler_vector;
  /**
   * A FiniteElement object which is only needed in 3D, since it knows how to
   * reorder shape functions/DoFs on non-standard faces. This is used to
   * reorder support points in the same way. We could make this a pointer to
   * prevent construction in 1D and 2D, but since memory and time requirements
   * are not particularly high this seems unnecessary at the moment.
   */
  SmartPointer<const FiniteElement<dim,spacedim>, MappingFEField<dim,spacedim,DH,VECTOR> > fe;


  /**
   * Pointer to the DoFHandler to which the mapping vector is associated.
   */
  SmartPointer<const DH,MappingFEField<dim,spacedim,DH,VECTOR> >euler_dof_handler;



private:
  /**
   * Update internal degrees of freedom.
   */
  void update_internal_dofs(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                            const typename MappingFEField<dim, spacedim>::InternalData &data) const;

  /**
   * Reimplemented from Mapping. See the documentation of the base class for
   * detailed information.
   */
  virtual void
  compute_shapes_virtual (const std::vector<Point<dim> > &unit_points,
                          typename MappingFEField<dim, spacedim>::InternalData &data) const;

  /**
   * Reimplemented from Mapping. See the documentation of the base class for
   * detailed information.
   */
  virtual UpdateFlags
  update_once (const UpdateFlags in) const;

  /**
   * Reimplemented from Mapping. See the documentation of the base class for
   * detailed information.
   */
  virtual UpdateFlags
  update_each (const UpdateFlags in) const;

  /**
   * Reimplemented from Mapping. See the documentation of the base class for
   * detailed information.
   */
  void
  compute_data (const UpdateFlags      update_flags,
                const Quadrature<dim>  &q,
                const unsigned int     n_original_q_points,
                InternalData           &data) const;

  /**
   * Reimplemented from Mapping. See the documentation of the base class for
   * detailed information.
   */
  void
  compute_face_data (const UpdateFlags      update_flags,
                     const Quadrature<dim>  &q,
                     const unsigned int     n_original_q_points,
                     InternalData           &data) const;

  /**
   * Reimplemented from Mapping. See the documentation of the base class for
   * detailed information.
   */
  virtual
  InternalData *
  get_data (const UpdateFlags,
            const Quadrature<dim> &quadrature) const;

  /**
   * Reimplemented from Mapping. See the documentation of the base class for
   * detailed information.
   */
  virtual
  typename Mapping<dim,spacedim>::InternalDataBase *
  get_face_data (const UpdateFlags flags,
                 const Quadrature<dim-1>& quadrature) const;

  /**
   * Reimplemented from Mapping. See the documentation of the base class for
   * detailed information.
   */
  virtual
  typename Mapping<dim,spacedim>::InternalDataBase *
  get_subface_data (const UpdateFlags flags,
                    const Quadrature<dim-1>& quadrature) const;


  /*
   * Which components to use for the mapping.
   */
  const ComponentMask fe_mask;


  /**
   * Mapping between indices in the FE space and the real space. This vector
   * contains one index for each component of the finite element space. If the
   * index is one for which the ComponentMask which is used to construct this
   * element is false, then numbers::invalid_unsigned_int is returned,
   * otherwise the component in real space is returned. For example, if we
   * construct the mapping using ComponentMask(spacedim, true), then this
   * vector contains {0,1,2} in spacedim = 3.
   */
  std::vector<unsigned int> fe_to_real;


  /**
   * Declare other MappingFEField classes friends.
   */
  template <int,int,class,class> friend class MappingFEField;
};

/*@}*/

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN


template<int dim, int spacedim, class DH, class VECTOR>
inline
const double &
MappingFEField<dim,spacedim,DH,VECTOR>::InternalData::shape (const unsigned int qpoint,
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
MappingFEField<dim,spacedim,DH,VECTOR>::InternalData::shape (const unsigned int qpoint,
    const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_values.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_values.size()));
  return shape_values [qpoint*n_shape_functions + shape_nr];
}


template<int dim, int spacedim, class DH, class VECTOR>
inline
const Tensor<1,dim> &
MappingFEField<dim,spacedim,DH,VECTOR>::InternalData::derivative (const unsigned int qpoint,
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
MappingFEField<dim,spacedim,DH,VECTOR>::InternalData::derivative (const unsigned int qpoint,
    const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_derivatives.size()));
  return shape_derivatives [qpoint*n_shape_functions + shape_nr];
}


template <int dim, int spacedim, class DH, class VECTOR>
inline
const Tensor<2,dim> &
MappingFEField<dim,spacedim,DH,VECTOR>::InternalData::second_derivative (const unsigned int qpoint,
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
MappingFEField<dim,spacedim,DH,VECTOR>::InternalData::second_derivative (const unsigned int qpoint,
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
MappingFEField<dim,spacedim,DH,VECTOR>::preserves_vertex_locations () const
{
  return false;
}




#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
