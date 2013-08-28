// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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

#ifndef __deal2__hp_dof_levels_h
#define __deal2__hp_dof_levels_h


#include <deal.II/base/config.h>
#include <deal.II/hp/dof_objects.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace hp
  {
    template <int dim>
    class DoFLevel;

    /**
     * Store the indices of the degrees of freedom that are located on
     * objects of dimension @p structdim.
     *
     * The things we store here are very similar to what is stored in the
     * internal::DoFHandler::DoFLevel class hierarchy (see there for more
     * information, in particular on the layout of the class hierarchy,
     * and the use of file names). There are two main
     * differences, discussed in the following subsections. In addition to
     * the data already stored by the internal::DoFHandler::DoFLevel
     * classes, we also have to store which finite element each cell
     * uses.
     *
     *
     * <h4>Offset computations</h4>
     *
     * For hp methods, not all cells may use the same finite element, and
     * it is consequently more complicated to determine where the DoF
     * indices for a given line, quad, or hex are stored. As described in
     * the documentation of the internal::DoFHandler::DoFLevel class, we
     * can compute the location of the first line DoF, for example, by
     * calculating the offset as <code>line_index *
     * dof_handler.get_fe().dofs_per_line</code>. This of course doesn't
     * work any more if different lines may have different numbers of
     * degrees of freedom associated with them. Consequently, rather than
     * using this simple multiplication, each of the lines.dofs, quads.dofs,
     * and hexes.dofs arrays has an associated array lines.dof_offsets,
     * quads.dof_offsets, and hexes.dof_offsets. The data corresponding to a
     * line then starts at index <code>line_dof_offsets[line_index]</code>
     * within the <code>line_dofs</code> array.
     *
     *
     * <h4>Multiple data sets per object</h4>
     *
     * If an object corresponds to a cell, the global dof indices of this
     * cell are stored at the location indicated above in sequential
     * order.
     *
     * However, if two adjacent cells use different finite elements, then
     * the face that they share needs to store DoF indices for both
     * involved finite elements. While faces therefore have to have at
     * most two sets of DoF indices, it is easy to see that vertices for
     * example can have as many sets of DoF indices associated with them
     * as there are adjacent cells, and the same holds for lines in 3d.
     *
     * Consequently, for objects that have a lower dimensionality than
     * cells, we have to store a map from the finite element index to the
     * set of DoF indices associated. Since real sets are typically very
     * inefficient to store, and since most of the time we expect the
     * number of individual keys to be small (frequently, adjacent cells
     * will have the same finite element, and only a single entry will
     * exist in the map), what we do is instead to store a linked list. In
     * this format, the first entry starting at position
     * <code>lines.dofs[lines.dof_offsets[line_index]]</code> will denote
     * the finite element index of the set of DoF indices following; after
     * this set, we will store the finite element index of the second set
     * followed by the corresponding DoF indices; and so on. Finally, when
     * all finite element indices adjacent to this object have been
     * covered, we write a -1 to indicate the end of the list.
     *
     * Access to this kind of data, as well as the distinction between
     * cells and objects of lower dimensionality are encoded in the
     * accessor functions, DoFObjects::set_dof_index() and
     * DoFLevel::get_dof_index() They are able to traverse this
     * list and pick out or set a DoF index given the finite element index
     * and its location within the set of DoFs corresponding to this
     * finite element.
     *
     *
     * @ingroup hp
     * @author Wolfgang Bangerth, 1998, 2006, Oliver Kayser-Herold 2003.
     */
    template <int dim>
    class DoFLevel
    {
    public:
      /**
       *  Indices specifying the finite
       *  element of hp::FECollection to use
       *  for the different cells on the current level. The
       *  meaning what a cell is, is
       *  dimension specific, therefore also
       *  the length of this vector depends
       *  on the dimension: in one dimension,
       *  the length of this vector equals
       *  the length of the @p lines vector,
       *  in two dimensions that of the @p
       *  quads vector, etc. The vector stores one element per cell
       *  since the actiev_fe_index is unique for cells.
       */
      std::vector<unsigned int> active_fe_indices;

      /**
       *  Store the dof-indices and related data of
       *  the cells on the current level corresponding to this object.
       */
      internal::hp::DoFObjects<dim> dof_object;

      /**
       * Determine an estimate for the
       * memory consumption (in bytes)
       * of this object.
       */
      std::size_t memory_consumption () const;
    };
  } // namespace hp

} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
