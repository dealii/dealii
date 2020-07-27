// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#ifndef dealii_non_local_dof_handler_h
#define dealii_non_local_dof_handler_h


#include <deal.II/base/config.h>

#include <deal.II/base/subscriptor.h>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim, bool>
class DoFCellAccessor;
#endif

/**
 * This class is used to enumerate non local degrees of freedom. Its default
 * implementation does nothing, since in general FiniteElement spaces only
 * define degrees of freedom on vertices, edges, faces, or cells. There are
 * cases, however, in which a FiniteElement space may define some non-local
 * basis functions which are non-zero on a given cell, even if the basis
 * function cannot be attributed to local subobjects of the given cell.
 *
 * In these cases, the DoFHandler class needs to know how to distribute these
 * extra degrees of freedom, and since it cannot do so on its own, it asks
 * the FiniteElement to return a NonLocalDoFHandler class that knows how to
 * enumerate and distribute these non local degrees of freedom.
 *
 * This class is intended as a base class for all FiniteElement spaces that
 * need to enumerate non local degrees of freedom.
 *
 * @ingroup dofs
 */
template <int dim, int spacedim = dim>
class NonLocalDoFHandler : public Subscriptor
{
public:
  /**
   * Make the dimension available in function templates.
   */
  static const unsigned int dimension = dim;

  /**
   * Make the space dimension available in function templates.
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * Standard constructor, not initializing any data.
   */
  NonLocalDoFHandler() = default;

  /**
   * Destructor.
   */
  virtual ~NonLocalDoFHandler() = default;

  /**
   * Copy operator. NonLocalDoFHandler objects may be large and expensive.
   * They should not be copied, in particular not by accident, but
   * rather deliberately constructed. As a consequence, this operator
   * is explicitly removed from the interface of this class.
   */
  NonLocalDoFHandler &
  operator=(const NonLocalDoFHandler &) = delete;

  /**
   * Return the non local dof indices associated to the current cell, for
   * active cell iterators.
   */
  virtual std::vector<types::global_dof_index>
  get_non_local_dof_indices(
    const DoFCellAccessor<dim, spacedim, false> &) const;

  /**
   * Return the non local dof indices associated to the current cell, for
   * level cell iterators.
   */
  virtual std::vector<types::global_dof_index>
  get_non_local_dof_indices(const DoFCellAccessor<dim, spacedim, true> &) const;

  /**
   * Return the number of non local dof indices that are required in addition to
   * the local ones.
   */
  virtual types::global_dof_index
  n_additional_non_local_dofs() const;

  /**
   * At any given moment, only one NonLocalDoFHandler can be used. Throw an
   * exception if two cells define an active FiniteElement for which
   * get_non_local_dof_handler() return a different object.
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcDifferentNonLocalDoFHandlers);
};



// ----------------------------------------------------------------------------
template <int dim, int spacedim>
inline std::vector<types::global_dof_index>
NonLocalDoFHandler<dim, spacedim>::get_non_local_dof_indices(
  const DoFCellAccessor<dim, spacedim, true> &) const
{
  // By default, we return an empty vector.
  return std::vector<types::global_dof_index>();
}



template <int dim, int spacedim>
inline std::vector<types::global_dof_index>
NonLocalDoFHandler<dim, spacedim>::get_non_local_dof_indices(
  const DoFCellAccessor<dim, spacedim, false> &) const
{
  // By default, we return an empty vector.
  return std::vector<types::global_dof_index>();
}



template <int dim, int spacedim>
inline types::global_dof_index
NonLocalDoFHandler<dim, spacedim>::n_additional_non_local_dofs() const
{
  return 0;
}

DEAL_II_NAMESPACE_CLOSE

#endif
