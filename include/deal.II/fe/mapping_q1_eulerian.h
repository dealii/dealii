// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_mapping_q1_eulerian_h
#define dealii_mapping_q1_eulerian_h

#include <deal.II/base/config.h>

#include <deal.II/base/observer_pointer.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/mapping_q.h>

#include <boost/container/small_vector.hpp>

#include <array>

DEAL_II_NAMESPACE_OPEN

template <typename>
class Vector;


/**
 * @addtogroup mapping
 * @{
 */

/**
 * This class provides a mapping that adds to the location of each cell
 * a $d$-linear displacement field. (The generalization to higher order
 * polynomials is provided in the MappingQEulerian class.) Each
 * cell is thus shifted in space by values given to the mapping through a
 * finite element field.
 *
 * <h3>Usage</h3>
 *
 * The two-argument constructor takes a reference to a displacement vector
 * and its DoFHandler. The finite element must have at least @p spacedim
 * components. The first @p spacedim components define the displacement;
 * a ComponentMask can select different components from a coupled solution.
 *
 * Each selected component must have exactly one primitive degree of freedom
 * at each vertex. Selected shape functions associated with other degrees of
 * freedom must vanish at vertices, since this mapping only reads vertex
 * coefficients. The other components can use different finite elements.
 * These requirements are met by @p spacedim continuous FE_Q components.
 *
 * The DoFHandler must have a finite element assigned before construction,
 * through DoFHandler::distribute_dofs(). The finite element must remain
 * unchanged for the lifetime of the mapping.
 *
 * An example is shown below:
 * @code
 *    FESystem<dim> fe(FE_Q<dim>(1), dim);
 *    DoFHandler<dim> flowfield_dof_handler(triangulation);
 *    flowfield_dof_handler.distribute_dofs(fe);
 *    Vector<double> displacement_field(flowfield_dof_handler.n_dofs());
 *    MappingQ1Eulerian<dim> mymapping(flowfield_dof_handler,
 *                                     displacement_field);
 * @endcode
 *
 * Note that since the vector of shift values and the dof handler are only
 * associated to this object at construction time, you have to make sure that
 * whenever you use this object, the given objects still represent valid data.
 *
 * To enable the use of the MappingQ1Eulerian class also in the context of
 * parallel codes using the PETSc or Trilinos wrapper classes, the type of
 * the vector can be specified as template parameter <tt>VectorType</tt>.
 *
 * For more information about the <tt>spacedim</tt> template parameter check
 * the documentation of FiniteElement or the one of Triangulation.
 */
template <int dim, typename VectorType = Vector<double>, int spacedim = dim>
class MappingQ1Eulerian : public MappingQ<dim, spacedim>
{
public:
  /**
   * Constructor.
   *
   * @param[in] euler_dof_handler A DoFHandler object that defines a finite
   * element space with at least @p spacedim components. The first @p spacedim
   * components are interpreted as displacements relative to the original
   * cell positions. Each selected component must satisfy the vertex
   * requirements described in the class documentation. The finite element
   * must already be assigned to the DoFHandler and remain unchanged for the
   * lifetime of the mapping.
   * @param[in] euler_vector A finite element function in the space defined by
   * the first argument. The first @p spacedim components define the displacement.
   */
  MappingQ1Eulerian(const DoFHandler<dim, spacedim> &euler_dof_handler,
                    const VectorType                &euler_vector);

  /**
   * Construct a mapping using the @p spacedim components selected by
   * @p component_mask, in increasing component order. Each selected component
   * must have exactly one primitive degree of freedom at each vertex. Other
   * components are ignored. Selected shape functions associated with DoFs
   * away from vertices must vanish at vertices.
   * An empty mask selects the first @p spacedim components, including in a
   * field with more components. The DoFHandler must have a finite element
   * assigned before construction, and that finite element must remain
   * unchanged for the lifetime of the mapping.
   */
  MappingQ1Eulerian(const DoFHandler<dim, spacedim> &euler_dof_handler,
                    const VectorType                &euler_vector,
                    const ComponentMask             &component_mask);


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
   * Always returns @p false because MappingQ1Eulerian does not in general
   * preserve vertex locations (unless the translation vector happens to
   * provide for zero displacements at vertex locations).
   */
  virtual bool
  preserves_vertex_locations() const override;

  /**
   * Exception.
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
   * Compute the support points of the mapping. For the current class, these
   * are the vertices, as obtained by calling Mapping::get_vertices(). See the
   * documentation of MappingQ::compute_mapping_support_points() for
   * more information.
   */
  virtual void
  compute_mapping_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    boost::container::small_vector<Point<spacedim>, 200>       &points)
    const override;

  /**
   * Reference to the vector of shifts.
   */
  ObserverPointer<const VectorType,
                  MappingQ1Eulerian<dim, VectorType, spacedim>>
    euler_transform_vectors;

  /**
   * Pointer to the DoFHandler to which the mapping vector is associated.
   */
  ObserverPointer<const DoFHandler<dim, spacedim>,
                  MappingQ1Eulerian<dim, VectorType, spacedim>>
    shiftmap_dof_handler;

  /**
   * Mask of the selected displacement components, with one entry per finite
   * element component and exactly @p spacedim entries selected.
   */
  const ComponentMask component_mask;

private:
  /**
   * Return the vertex DoF indices of the selected displacement components.
   */
  std::array<unsigned int, spacedim>
  get_vertex_dof_indices() const;
};

/** @} */

/*----------------------------------------------------------------------*/

#ifndef DOXYGEN

template <int dim, typename VectorType, int spacedim>
inline bool
MappingQ1Eulerian<dim, VectorType, spacedim>::preserves_vertex_locations() const
{
  return false;
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
