// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2016 by the deal.II authors
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

#ifndef dealii__mapping_q1_eulerian_h
#define dealii__mapping_q1_eulerian_h

#include <deal.II/base/config.h>
#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/dofs/dof_handler.h>

DEAL_II_NAMESPACE_OPEN

template <typename> class Vector;


/*!@addtogroup mapping */
/*@{*/

/**
 * Eulerian mapping of general unit cells by $d$-linear shape functions. Each
 * cell is thus shifted in space by values given to the mapping through a
 * finite element field.
 *
 * <h3>Usage</h3>
 *
 * The constructor of this class takes two arguments: a reference to the
 * vector that defines the mapping from the reference configuration to the
 * current configuration and a reference to the DoFHandler. The vector should
 * then represent a (flattened out version of a) vector valued field defined
 * at nodes defined by the the DoFHandler, where the number of components of
 * the vector field equals the number of space dimensions. Thus, the
 * DoFHandler shall operate on a finite element that has as many components as
 * space dimensions. As an additional requirement, we impose that it have as
 * many degree of freedom per vertex as there are space dimensions; since this
 * object only evaluates the finite element field at the vertices, the values
 * of all other degrees of freedom (not associated to vertices) are ignored.
 * These requirements are met if the finite element which the given DoFHandler
 * operates on is constructed as a system element (FESystem) from @p dim
 * continuous FE_Q() objects.
 *
 * In many cases, the shift vector will also be the solution vector of the
 * problem under investigation. If this is not the case (i.e. the number of
 * components of the solution variable is not equal to the space dimension,
 * e.g. for scalar problems in <tt>dim>1</tt> where the Eulerian coordinates
 * only give a background field) or for coupled problems where more variables
 * are computed than just the flow field), then a different DoFHandler has to
 * be set up on the given triangulation, and the shift vector has then to be
 * associated to it.
 *
 * An example is shown below:
 * @code
 *    FESystem<dim> fe(FE_Q<dim>(1), dim);
 *    DoFHandler<dim> flowfield_dof_handler(triangulation);
 *    flowfield_dof_handler.distribute_dofs(fe);
 *    Vector<double> map_points(flowfield_dof_handler.n_dofs());
 *    MappingQ1Eulerian<dim> mymapping(map_points, flowfield_dof_handler);
 * @endcode
 *
 * Note that since the vector of shift values and the dof handler are only
 * associated to this object at construction time, you have to make sure that
 * whenever you use this object, the given objects still represent valid data.
 *
 * To enable the use of the MappingQ1Eulerian class also in the context of
 * parallel codes using the PETSc wrapper classes, the type of the vector can
 * be specified as template parameter <tt>EulerVectorType</tt> Not specifying
 * this template argument in applications using the PETSc vector classes leads
 * to the construction of a copy of the vector which is not acccessible
 * afterwards!
 *
 * For more information about the <tt>spacedim</tt> template parameter check
 * the documentation of FiniteElement or the one of Triangulation.
 *
 * @author Michael Stadler, 2001
 */
template <int dim, typename VectorType = Vector<double>, int spacedim=dim >
class MappingQ1Eulerian : public MappingQGeneric<dim,spacedim>
{
public:

  /**
   * Constructor. It takes a <tt>Vector<double> &</tt> as its first argument
   * to specify the transformation of the whole problem from the reference to
   * the current configuration. The organization of the elements in the @p
   * Vector must follow the concept how deal.II stores solutions that are
   * associated to a triangulation.  This is automatically the case if the @p
   * Vector represents the solution of the previous step of a nonlinear
   * problem. Alternatively, the @p Vector can be initialized by
   * <tt>DoFAccessor::set_dof_values()</tt>.
   */
  MappingQ1Eulerian (const VectorType  &euler_transform_vectors,
                     const DoFHandler<dim,spacedim> &shiftmap_dof_handler);

  /**
   * Return the mapped vertices of the cell. For the current class, this
   * function does not use the support points from the geometry of the current
   * cell but instead evaluates an externally given displacement field in
   * addition to the geometry of the cell.
   */
  virtual
  std_cxx11::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
  get_vertices (const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;

  /**
   * Return a pointer to a copy of the present object. The caller of this copy
   * then assumes ownership of it.
   */
  virtual
  MappingQ1Eulerian<dim,VectorType,spacedim> *clone () const;

  /**
   * Always returns @p false because MappingQ1Eulerian does not in general
   * preserve vertex locations (unless the translation vector happens to
   * provide for zero displacements at vertex locations).
   */
  bool preserves_vertex_locations () const;

  /**
   * Exception.
   */
  DeclException0 (ExcInactiveCell);



protected:
  /**
   * Compute mapping-related information for a cell. See the documentation of
   * Mapping::fill_fe_values() for a discussion of purpose, arguments, and
   * return value of this function.
   *
   * This function overrides the function in the base class since we cannot
   * use any cell similarity for this class.
   */
  virtual
  CellSimilarity::Similarity
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                  const CellSimilarity::Similarity                           cell_similarity,
                  const Quadrature<dim>                                     &quadrature,
                  const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                  internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const;

  /**
   * Compute the support points of the mapping. For the current class, these
   * are the vertices, as obtained by calling Mapping::get_vertices(). See the
   * documentation of MappingQGeneric::compute_mapping_support_points() for
   * more information.
   */
  virtual
  std::vector<Point<spacedim> >
  compute_mapping_support_points(const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;

  /**
   * Reference to the vector of shifts.
   */
  SmartPointer<const VectorType, MappingQ1Eulerian<dim,VectorType,spacedim> > euler_transform_vectors;

  /**
   * Pointer to the DoFHandler to which the mapping vector is associated.
   */
  SmartPointer<const DoFHandler<dim,spacedim>,MappingQ1Eulerian<dim,VectorType,spacedim> > shiftmap_dof_handler;
};

/*@}*/

/*----------------------------------------------------------------------*/

#ifndef DOXYGEN

template <int dim, typename VectorType, int spacedim>
inline
bool
MappingQ1Eulerian<dim,VectorType,spacedim>::preserves_vertex_locations () const
{
  return false;
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
