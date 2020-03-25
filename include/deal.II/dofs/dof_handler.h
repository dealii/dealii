// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_dof_handler_h
#define dealii_dof_handler_h



#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/block_info.h>
#include <deal.II/dofs/deprecated_function_map.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_faces.h>
#include <deal.II/dofs/dof_handler_base.h>
#include <deal.II/dofs/dof_iterator_selector.h>
#include <deal.II/dofs/dof_levels.h>
#include <deal.II/dofs/number_cache.h>

#include <deal.II/hp/fe_collection.h>

#include <boost/serialization/split_member.hpp>

#include <map>
#include <memory>
#include <set>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * Given a triangulation and a description of a finite element, this
 * class enumerates degrees of freedom on all vertices, edges, faces,
 * and cells of the triangulation. As a result, it also provides a
 * <i>basis</i> for a discrete space $V_h$ whose elements are finite
 * element functions defined on each cell by a FiniteElement object.
 * This class satisfies the
 * @ref ConceptMeshType "MeshType concept"
 * requirements.
 *
 * It is first used in the step-2 tutorial program.
 *
 * For each vertex, line, quad, etc, this class stores a list of the indices
 * of degrees of freedom living on this object. These indices refer to the
 * unconstrained degrees of freedom, i.e. constrained degrees of freedom are
 * numbered in the same way as unconstrained ones, and are only later
 * eliminated.  This leads to the fact that indices in global vectors and
 * matrices also refer to all degrees of freedom and some kind of condensation
 * is needed to restrict the systems of equations to the unconstrained degrees
 * of freedom only. The actual layout of storage of the indices is described
 * in the dealii::internal::DoFHandlerImplementation::DoFLevel class
 * documentation.
 *
 * The class offers iterators to traverse all cells, in much the same way as
 * the Triangulation class does. Using the begin() and end() functions (and
 * companions, like begin_active()), one can obtain iterators to walk over
 * cells, and query the degree of freedom structures as well as the
 * triangulation data. These iterators are built on top of those of the
 * Triangulation class, but offer the additional information on degrees of
 * freedom functionality compared to pure triangulation iterators. The order
 * in which dof iterators are presented by the <tt>++</tt> and <tt>\--</tt>
 * operators is the same as that for the corresponding iterators traversing
 * the triangulation on which this DoFHandler is constructed.
 *
 * The <tt>spacedim</tt> parameter has to be used if one wants to solve
 * problems on surfaces. If not specified, this parameter takes the default
 * value <tt>=dim</tt> implying that we want to solve problems in a domain
 * whose dimension equals the dimension of the space in which it is embedded.
 *
 *
 * <h3>Distribution of indices for degrees of freedom</h3>
 *
 * The degrees of freedom (`dofs') are distributed on the given triangulation
 * by the function distribute_dofs(). It gets passed a finite element object
 * describing how many degrees of freedom are located on vertices, lines, etc.
 * It traverses the triangulation cell by cell and numbers the dofs of that
 * cell if not yet numbered. For non-multigrid algorithms, only active cells
 * are considered. Active cells are defined to be those cells which have no
 * children, i.e. they are the most refined ones.
 *
 * Since the triangulation is traversed starting with the cells of the
 * coarsest active level and going to more refined levels, the lowest numbers
 * for dofs are given to the largest cells as well as their bounding lines and
 * vertices, with the dofs of more refined cells getting higher numbers.
 *
 * This numbering implies very large bandwidths of the resulting matrices and
 * is thus vastly suboptimal for some solution algorithms. For this reason,
 * the DoFRenumbering class offers several algorithms to reorder the dof
 * numbering according. See there for a discussion of the implemented
 * algorithms.
 *
 *
 * <h3>Interaction with distributed meshes</h3>
 *
 * Upon construction, this class takes a reference to a triangulation object.
 * In most cases, this will be a reference to an object of type Triangulation,
 * i.e. the class that represents triangulations that entirely reside on a
 * single processor. However, it can also be of type
 * parallel::distributed::Triangulation (see, for example, step-32, step-40
 * and in particular the
 * @ref distributed
 * module) in which case the DoFHandler object will proceed to only manage
 * degrees of freedom on locally owned and ghost cells. This process is
 * entirely transparent to the used.
 *
 *
 * <h3>User defined renumbering schemes</h3>
 *
 * The DoFRenumbering class offers a number of renumbering schemes like the
 * Cuthill-McKee scheme. Basically, the function sets up an array in which for
 * each degree of freedom we store the new index this DoF should have after
 * renumbering. Using this array, the renumber_dofs() function of the present
 * class is called, which actually performs the change from old DoF indices to
 * the ones given in the array. In some cases, however, a user may want to
 * compute their own renumbering order; in this case, one can allocate an array
 * with one element per degree of freedom and fill it with the number that the
 * respective degree of freedom shall be assigned. This number may, for
 * example, be obtained by sorting the support points of the degrees of
 * freedom in downwind direction.  Then call the
 * <tt>renumber_dofs(vector<types::global_dof_index>)</tt> function with the
 * array, which converts old into new degree of freedom indices.
 *
 *
 * <h3>Serializing (loading or storing) DoFHandler objects</h3>
 *
 * Like many other classes in deal.II, the DoFHandler class can stream its
 * contents to an archive using BOOST's serialization facilities. The data so
 * stored can later be retrieved again from the archive to restore the
 * contents of this object. This facility is frequently used to save the state
 * of a program to disk for possible later resurrection, often in the context
 * of checkpoint/restart strategies for long running computations or on
 * computers that aren't very reliable (e.g. on very large clusters where
 * individual nodes occasionally fail and then bring down an entire MPI job).
 *
 * The model for doing so is similar for the DoFHandler class as it is for the
 * Triangulation class (see the section in the general documentation of that
 * class). In particular, the load() function does not exactly restore the
 * same state as was stored previously using the save() function. Rather, the
 * function assumes that you load data into a DoFHandler object that is
 * already associated with a triangulation that has a content that matches the
 * one that was used when the data was saved. Likewise, the load() function
 * assumes that the current object is already associated with a finite element
 * object that matches the one that was associated with it when data was
 * saved; the latter can be achieved by calling DoFHandler::distribute_dofs()
 * using the same kind of finite element before re-loading data from the
 * serialization archive.
 *
 * @ingroup dofs
 * @author Wolfgang Bangerth, Markus Buerg, Timo Heister, Guido Kanschat
 * @date 1998, 1999, 2000, 2012
 *
 * @note Task is delegated to the base class DoFHandlerBase.
 */
template <int dim, int spacedim = dim>
class DoFHandler
  : public DoFHandlerBase<dim, spacedim, DoFHandler<dim, spacedim>>
{
public:
  /**
   * Make the type of this DoFHandler available in function templates.
   */
  static const bool is_hp_dof_handler = false;

  /**
   * Standard constructor, not initializing any data. After constructing an
   * object with this constructor, use initialize() to make a valid
   * DoFHandler.
   */
  DoFHandler();

  /**
   * Constructor. Take @p tria as the triangulation to work on.
   */
  DoFHandler(const Triangulation<dim, spacedim> &tria);

  /**
   * Copy constructor. DoFHandler objects are large and expensive.
   * They should not be copied, in particular not by accident, but
   * rather deliberately constructed. As a consequence, this constructor
   * is explicitly removed from the interface of this class.
   */
  DoFHandler(const DoFHandler &) = delete;

  /**
   * Copy operator. DoFHandler objects are large and expensive.
   * They should not be copied, in particular not by accident, but
   * rather deliberately constructed. As a consequence, this operator
   * is explicitly removed from the interface of this class.
   */
  DoFHandler &
  operator=(const DoFHandler &) = delete;
};


DEAL_II_NAMESPACE_CLOSE

#endif
