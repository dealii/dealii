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


#ifndef dealii__mapping_q_eulerian_h
#define dealii__mapping_q_eulerian_h

#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>


DEAL_II_NAMESPACE_OPEN

template <typename> class Vector;


/*!@addtogroup mapping */
/*@{*/

/**
 * This class is an extension of the MappingQ1Eulerian class to higher order
 * Qp mappings.  It is useful when one wants to calculate shape function
 * information on a domain that is deforming as the computation proceeds.
 *
 * <h3>Usage</h3>
 *
 * The constructor of this class takes three arguments: the polynomial degree
 * of the desire Qp mapping, a reference to the vector that defines the
 * mapping from the initial configuration to the current configuration, and a
 * reference to the DoFHandler. The most common case is to use the solution
 * vector for the problem under consideration as the shift vector. The key
 * requirement is that the number of components of the given vector field be
 * equal to (or possibly greater than) the number of space dimensions. If
 * there are more components than space dimensions (for example, if one is
 * working with a coupled problem where there are additional solution
 * variables), the first <tt>dim</tt> components are assumed to represent the
 * displacement field, and the remaining components are ignored.  If this
 * assumption does not hold one may need to set up a separate DoFHandler on
 * the triangulation and associate the desired shift vector to it.
 *
 * Typically, the DoFHandler operates on a finite element that is constructed
 * as a system element (FESystem) from continuous FE_Q() objects. An example
 * is shown below:
 * @code
 *    FESystem<dim> fe(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1);
 *    DoFHandler<dim> dof_handler(triangulation);
 *    dof_handler.distribute_dofs(fe);
 *    Vector<double> displacement_field(dof_handler.n_dofs());
 *    // ... compute displacement field somehow...
 *    MappingQEulerian<dim> q2_mapping(2, dof_handler, displacement_field);
 * @endcode
 *
 * In this example, our element consists of <tt>(dim+1)</tt> components. Only
 * the first <tt>dim</tt> components will be used, however, to define the Q2
 * mapping.  The remaining components are ignored.
 *
 * Note that it is essential to call the distribute_dofs(...) function before
 * constructing a mapping object.
 *
 * Also note that since the vector of shift values and the dof handler are
 * only associated to this object at construction time, you have to make sure
 * that whenever you use this object, the given objects still represent valid
 * data.
 *
 * To enable the use of the MappingQ1Eulerian class also in the context of
 * parallel codes using the PETSc wrapper classes, the type of the vector can
 * be specified as template parameter <tt>EulerVectorType</tt> Not specifying
 * this template argument in applications using the PETSc vector classes leads
 * to the construction of a copy of the vector which is not accessible
 * afterwards!
 *
 * @author Joshua White, 2008
 */
template <int dim, typename VectorType = Vector<double>, int spacedim=dim >
class MappingQEulerian : public MappingQ<dim, spacedim>
{
public:
  /**
   * Constructor.
   *
   * @param[in] degree The polynomial degree of the desired $Q_p$ mapping.
   * @param[in] euler_dof_handler A DoFHandler object that defines a finite
   * element space. This space needs to have at least dim components and the
   * first dim components of the space will be considered displacements
   * relative to the original positions of the cells of the triangulation.
   * @param[in] euler_vector A finite element function in the space defined by
   * the second argument. The first dim components of this function will be
   * interpreted as the displacement we use in defining the mapping, relative
   * to the location of cells of the underlying triangulation.
   */
  MappingQEulerian (const unsigned int             degree,
                    const DoFHandler<dim,spacedim> &euler_dof_handler,
                    const VectorType               &euler_vector);

  /**
   * @deprecated Use the constructor with the reverse order of second and
   * third argument.
   */
  MappingQEulerian (const unsigned int             degree,
                    const VectorType               &euler_vector,
                    const DoFHandler<dim,spacedim> &euler_dof_handler) DEAL_II_DEPRECATED;

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
  Mapping<dim,spacedim> *clone () const;

  /**
   * Always returns @p false because MappingQ1Eulerian does not in general
   * preserve vertex locations (unless the translation vector happens to
   * provide for zero displacements at vertex locations).
   */
  bool preserves_vertex_locations () const;

  /**
   * Exception
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
   * Reference to the vector of shifts.
   */
  SmartPointer<const VectorType, MappingQEulerian<dim,VectorType,spacedim> > euler_vector;

  /**
   * Pointer to the DoFHandler to which the mapping vector is associated.
   */
  SmartPointer<const DoFHandler<dim,spacedim>,MappingQEulerian<dim,VectorType,spacedim> > euler_dof_handler;


private:

  /**
   * A class derived from MappingQGeneric that provides the generic mapping
   * with support points on boundary objects so that the corresponding Q3
   * mapping ends up being C1.
   */
  class MappingQEulerianGeneric : public MappingQGeneric<dim,spacedim>
  {
  public:

    /**
     * Constructor.
     */
    MappingQEulerianGeneric (const unsigned int degree,
                             const MappingQEulerian<dim,VectorType,spacedim> &mapping_q_eulerian);

    /**
     * Return the mapped vertices of the cell. For the current class, this
     * function does not use the support points from the geometry of the
     * current cell but instead evaluates an externally given displacement
     * field in addition to the geometry of the cell.
     */
    virtual
    std_cxx11::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
    get_vertices (const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;

    /**
     * Compute the positions of the support points in the current
     * configuration. See the documentation of
     * MappingQGeneric::compute_mapping_support_points() for more information.
     */
    virtual
    std::vector<Point<spacedim> >
    compute_mapping_support_points(const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;

  private:
    /**
     * Reference to the surrounding object off of which we live.
     */
    const MappingQEulerian<dim,VectorType,spacedim> &mapping_q_eulerian;


    /**
     * Special quadrature rule used to define the support points in the
     * reference configuration.
     */

    class SupportQuadrature : public Quadrature<dim>
    {
    public:
      /**
       * Constructor, with an argument defining the desired polynomial degree.
       */

      SupportQuadrature (const unsigned int map_degree);

    };

    /**
     * A member variable holding the quadrature points in the right order.
     */
    const SupportQuadrature support_quadrature;

    /**
     * FEValues object used to query the the given finite element field at the
     * support points in the reference configuration.
     *
     * The variable is marked as mutable since we have to call
     * FEValues::reinit from compute_mapping_support_points, a function that
     * is 'const'.
     */
    mutable FEValues<dim,spacedim> fe_values;

    /**
     * A variable to guard access to the fe_values variable.
     */
    mutable Threads::Mutex fe_values_mutex;
  };

};

/*@}*/


/*----------------------------------------------------------------------*/

#ifndef DOXYGEN

template <int dim, typename VectorType, int spacedim>
inline
bool
MappingQEulerian<dim,VectorType,spacedim>::preserves_vertex_locations () const
{
  return false;
}

#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE


#endif // dealii__mapping_q_eulerian_h
