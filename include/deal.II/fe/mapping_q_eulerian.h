// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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


#ifndef __deal2__mapping_q_eulerian_h
#define __deal2__mapping_q_eulerian_h

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>


DEAL_II_NAMESPACE_OPEN


/*!@addtogroup mapping */
/*@{*/

/**
 * This class is an extension of the MappingQ1Eulerian
 * class to higher order Qp mappings.  It is useful
 * when one wants to calculate shape function information on
 * a domain that is deforming as the computation proceeds.
 *
 * <h3>Usage</h3>
 *
 * The constructor of this class takes three arguments: the polynomial
 * degree of the desire Qp mapping, a reference to
 * the vector that defines the mapping from the initial
 * configuration to the current configuration, and a reference to the
 * DoFHandler. The most common case is to use the solution
 * vector for the problem under consideration as the shift vector.
 * The key reqirement is that the number of components
 * of the given vector field be equal to (or possibly greater than) the
 * number of space dimensions. If there are more components than space
 * dimensions (for example, if one is working with a coupled problem
 * where there are additional solution variables), the
 * first <tt>dim</tt> components are assumed to represent the displacement
 * field, and the remaining components are ignored.  If this assumption
 * does not hold one may need to set up a separate DoFHandler on
 * the triangulation and associate the desired shift vector to it.
 *
 * Typically, the DoFHandler operates on a finite element that
 * is constructed as a system element (FESystem) from continuous FE_Q()
 * objects. An example is shown below:
 * @code
 *    FESystem<dim> fe(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1);
 *    DoFHandler<dim> dof_handler(triangulation);
 *    dof_handler.distribute_dofs(fe);
 *    Vector<double> soln_vector(dof_handler.n_dofs());
 *    MappingQEulerian<dim> q2_mapping(2,soln_vector,dof_handler);
 * @endcode
 *
 * In this example, our element consists of <tt>(dim+1)</tt> components.
 * Only the first <tt>dim</tt> components will be used, however, to define
 * the Q2 mapping.  The remaining components are ignored.
 *
 * Note that it is essential to call the distribute_dofs(...) function
 * before constructing a mapping object.
 *
 * Also note that since the vector of shift values and the dof handler are
 * only associated to this object at construction time, you have to
 * make sure that whenever you use this object, the given objects
 * still represent valid data.
 *
 * To enable the use of the MappingQ1Eulerian class also in the context
 * of parallel codes using the PETSc wrapper classes, the type of
 * the vector can be specified as template parameter <tt>EulerVectorType</tt>
 * Not specifying this template argument in applications using the PETSc
 * vector classes leads to the construction of a copy of the vector
 * which is not acccessible afterwards!
 *
 * @author Joshua White, 2008
 */
template <int dim, class VECTOR = Vector<double>, int spacedim=dim >
class MappingQEulerian : public MappingQ<dim, spacedim>
{
public:
  /**
   * Constructor. The first argument is the polynomical degree of the desired
   * Qp mapping.  It then takes a <tt>Vector<double> &</tt> to specify the
   * transformation of the domain from the reference to the current
   * configuration.  The organization of the elements in the @p Vector must
   * follow the concept how deal.II stores solutions that are associated to a
   * triangulation.  This is automatically the case if the @p Vector
   * represents the solution of the previous step of a nonlinear problem.
   * Alternatively, the @p Vector can be initialized by
   * <tt>DoFAccessor::set_dof_values()</tt>.
   */

  MappingQEulerian (const unsigned int     degree,
                    const VECTOR  &euler_vector,
                    const DoFHandler<dim,spacedim>  &euler_dof_handler);

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
   * Implementation of the interface in MappingQ. Overrides the function in
   * the base class, since we cannot use any cell similarity for this class.
   */
  virtual void
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                  const Quadrature<dim>                                     &quadrature,
                  typename Mapping<dim,spacedim>::InternalDataBase          &mapping_data,
                  typename std::vector<Point<spacedim> >                    &quadrature_points,
                  std::vector<double>                                       &JxW_values,
                  std::vector<DerivativeForm<1,dim,spacedim> >      &jacobians,
                  std::vector<DerivativeForm<2,dim,spacedim> >      &jacobian_grads,
                  std::vector<DerivativeForm<1,spacedim,dim> >      &inverse_jacobians,
                  std::vector<Point<spacedim> >                             &cell_normal_vectors,
                  CellSimilarity::Similarity                           &cell_similarity) const;

  /**
   * Reference to the vector of shifts.
   */

  SmartPointer<const VECTOR, MappingQEulerian<dim,VECTOR,spacedim> > euler_vector;

  /**
   * Pointer to the DoFHandler to which the mapping vector is associated.
   */
  SmartPointer<const DoFHandler<dim,spacedim>,MappingQEulerian<dim,VECTOR,spacedim> > euler_dof_handler;


private:

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
   * The variable is marked as mutable since we have to call FEValues::reinit
   * from compute_mapping_support_points, a function that is 'const'.
   */
  mutable FEValues<dim,spacedim> fe_values;

  /**
   * A variable to guard access to the fe_values variable.
   */
  mutable Threads::Mutex fe_values_mutex;

  /**
   * Compute the positions of the support points in the current configuration
   */
  virtual void compute_mapping_support_points(
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    std::vector<Point<spacedim> > &a) const;

};

/*@}*/


/*----------------------------------------------------------------------*/

#ifndef DOXYGEN

template <int dim, class VECTOR, int spacedim>
inline
bool
MappingQEulerian<dim,VECTOR,spacedim>::preserves_vertex_locations () const
{
  return false;
}

#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE


#endif // __deal2__mapping_q_eulerian_h

