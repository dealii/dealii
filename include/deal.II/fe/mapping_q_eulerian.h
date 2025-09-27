// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_mapping_q_eulerian_h
#define dealii_mapping_q_eulerian_h

#include <deal.II/base/config.h>

#include <deal.II/base/mutex.h>
#include <deal.II/base/observer_pointer.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/tria_iterator.h>

#include <boost/container/small_vector.hpp>


DEAL_II_NAMESPACE_OPEN

template <typename>
class Vector;


/**
 * @addtogroup mapping
 * @{
 */

/**
 * This class is an extension of the MappingQ1Eulerian class to higher order
 * $Q_p$ mappings.  It is useful when one wants to calculate shape function
 * information on a domain that is deforming as the computation proceeds.
 *
 * <h3>Usage</h3>
 *
 * The constructor of this class takes three arguments: the polynomial degree
 * of the desired Qp mapping, a reference to the vector that defines the
 * mapping from the initial configuration to the current configuration, and a
 * reference to the DoFHandler. The most common case is to use the solution
 * vector for the problem under consideration as the shift vector. The key
 * requirement is that the number of components of the given vector field must
 * be equal to (or possibly greater than) the number of space dimensions. If
 * there are more components than space dimensions (for example, if one is
 * working with a coupled problem where there are additional solution
 * variables), the first <tt>dim</tt> components are assumed to represent the
 * displacement field, and the remaining components are ignored.  If this
 * assumption does not hold one may need to set up a separate DoFHandler on
 * the triangulation and associate the desired shift vector to it.
 *
 * Typically, the DoFHandler operates on a finite element that is constructed
 * as a system element (FESystem) from continuous FE_Q objects. An example
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
 * To enable the use of the MappingQEulerian class also in the context of
 * parallel codes using the PETSc or Trilinos wrapper classes, the type
 * of the vector can be specified as template parameter <tt>VectorType</tt>.
 */
template <int dim, typename VectorType = Vector<double>, int spacedim = dim>
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
   * @param[in] level The multi-grid level at which the mapping will
   * be used. It is mainly used to check if the size of the @p euler_vector
   * is consistent with the @p euler_dof_handler .
   */
  MappingQEulerian(const unsigned int               degree,
                   const DoFHandler<dim, spacedim> &euler_dof_handler,
                   const VectorType                &euler_vector,
                   const unsigned int level = numbers::invalid_unsigned_int);

  /**
   * Return the mapped vertices of the cell. For the current class, this
   * function does not use the support points from the geometry of the current
   * cell but instead evaluates an externally given displacement field in
   * addition to the geometry of the cell.
   */
  virtual boost::container::small_vector<Point<spacedim>,
#ifndef _MSC_VER
                                         ReferenceCells::max_n_vertices<dim>()
#else
                                         GeometryInfo<dim>::vertices_per_cell
#endif
                                         >
  get_vertices(const typename Triangulation<dim, spacedim>::cell_iterator &cell)
    const override;

  /**
   * Return a pointer to a copy of the present object. The caller of this copy
   * then assumes ownership of it.
   */
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

  /**
   * Always return @p false because MappingQEulerian does not in general
   * preserve vertex locations (unless the translation vector happens to
   * provide zero displacements at vertex locations).
   */
  virtual bool
  preserves_vertex_locations() const override;

  // for documentation, see the Mapping base class
  virtual std::vector<Point<spacedim>>
  compute_mapping_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell)
    const override;

  /**
   * Exception which is thrown when the mapping is being evaluated at
   * non-active cell.
   */
  DeclException0(ExcInactiveCell);

protected:
  /**
   * Compute mapping-related information for a cell. See the documentation of
   * Mapping::fill_fe_values() for a discussion of purpose, arguments, and
   * return value of this function.
   *
   * This function overrides the function in the base class since we cannot
   * use any cell similarity for this class.
   */
  virtual CellSimilarity::Similarity
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim>                                      &quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  /**
   * Reference to the vector of shifts.
   */
  ObserverPointer<const VectorType, MappingQEulerian<dim, VectorType, spacedim>>
    euler_vector;

  /**
   * Pointer to the DoFHandler to which the mapping vector is associated.
   */
  ObserverPointer<const DoFHandler<dim, spacedim>,
                  MappingQEulerian<dim, VectorType, spacedim>>
    euler_dof_handler;

private:
  /**
   * Multigrid level at which the mapping is to be used.
   */
  const unsigned int level;

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
    SupportQuadrature(const unsigned int map_degree);
  };

  /**
   * A member variable holding the quadrature points in the right order.
   */
  const SupportQuadrature support_quadrature;

  /**
   * A MappingQ object, which is used by the fe_values
   * member variable to compute the undeformed mapping support
   * points, before adding any deformation.
   */
  const MappingQ<dim, spacedim> mapping_q;

  /**
   * FEValues object used to query the given finite element field at the
   * support points in the reference configuration.
   *
   * The variable is marked as mutable since we have to call
   * FEValues::reinit from compute_mapping_support_points, a function that
   * is 'const'.
   */
  mutable FEValues<dim, spacedim> fe_values;

  /**
   * A variable to guard access to the fe_values variable.
   */
  mutable Threads::Mutex fe_values_mutex;
};

/** @} */


/*----------------------------------------------------------------------*/

#ifndef DOXYGEN

template <int dim, typename VectorType, int spacedim>
inline bool
MappingQEulerian<dim, VectorType, spacedim>::preserves_vertex_locations() const
{
  return false;
}

#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE


#endif // dealii_mapping_q_eulerian_h
