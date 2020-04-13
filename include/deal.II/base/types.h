// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2018 by the deal.II authors
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

#ifndef dealii_types_h
#define dealii_types_h


#include <deal.II/base/config.h>

#include <cstdint>


DEAL_II_NAMESPACE_OPEN

/**
 * A namespace in which we declare alias for types used in deal.II, as well
 * as special values for these types.
 */
namespace types
{
  /**
   * The type used to denote subdomain_ids of cells.
   *
   * See the
   * @ref GlossSubdomainId "glossary"
   * for more information.
   *
   * There is a special value, numbers::invalid_subdomain_id that is used to
   * indicate an invalid value of this type.
   */
  using subdomain_id = unsigned int;

  /**
   * The type used for global indices of vertices.
   */
  using global_vertex_index = uint64_t;

  /**
   * An identifier that denotes the MPI type associated with
   * types::global_vertex_index.
   */
#define DEAL_II_VERTEX_INDEX_MPI_TYPE MPI_UINT64_T

  /**
   * The type used to denote the global index of degrees of freedom. This
   * type is then also used for querying the global *number* of degrees
   * of freedom, since the number is simply the largest index plus one.
   *
   * While in sequential computations the 4 billion indices of 32-bit unsigned
   * integers is plenty, parallel computations using (for example) the
   * parallel::distributed::Triangulation class can overflow this number and
   * consequently, deal.II chooses a larger integer type when
   * configured to use 64-bit indices.
   *
   * The data type always corresponds to an unsigned integer type.
   *
   * See the
   * @ref GlobalDoFIndex
   * page for guidance on when this type should or should not be used.
   */
#ifdef DEAL_II_WITH_64BIT_INDICES
  using global_dof_index = uint64_t;
#else
  using global_dof_index  = unsigned int;
#endif

  /**
   * An identifier that denotes the MPI type associated with
   * types::global_dof_index.
   */
#ifdef DEAL_II_WITH_64BIT_INDICES
#  define DEAL_II_DOF_INDEX_MPI_TYPE MPI_UINT64_T
#else
#  define DEAL_II_DOF_INDEX_MPI_TYPE MPI_UNSIGNED
#endif

  /**
   * The type used to denote the global index of a cell. This type
   * is then also used for querying the global *number* of cells in
   * a triangulation since the number is simply the largest index plus one.
   *
   * While in sequential computations the 4 billion indices of 32-bit unsigned
   * integers is plenty, parallel computations using (for example) the
   * parallel::distributed::Triangulation class can overflow this number and
   * consequently, deal.II chooses a larger integer type when
   * configured to use 64-bit indices.
   *
   * The data type always corresponds to an unsigned integer type.
   */
#ifdef DEAL_II_WITH_64BIT_INDICES
  using global_cell_index = uint64_t;
#else
  using global_cell_index = unsigned int;
#endif

  /**
   * The type used for coarse-cell ids. See the glossary
   * entry on @ref GlossCoarseCellId "coarse cell IDs" for more information.
   */
  using coarse_cell_id = global_cell_index;

  /**
   * The type used to denote boundary indicators associated with every piece
   * of the boundary and, in the case of meshes that describe manifolds in
   * higher dimensions, associated with every cell.
   *
   * There is a special value, numbers::internal_face_boundary_id that is used
   * to indicate an invalid value of this type and that is used as the
   * boundary indicator for faces that are in the interior of the domain and
   * therefore not part of any addressable boundary component.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  using boundary_id = unsigned int;

  /**
   * The type used to denote manifold indicators associated with every object
   * of the mesh.
   *
   * There is a special value, numbers::flat_manifold_id that is used to
   * indicate the standard cartesian manifold.
   *
   * @see
   * @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  using manifold_id = unsigned int;

  /**
   * The type used to denote material indicators associated with every cell.
   *
   * There is a special value, numbers::invalid_material_id that is used to
   * indicate an invalid value of this type.
   *
   * @see
   * @ref GlossMaterialId "Glossary entry on material indicators"
   */
  using material_id = unsigned int;
} // namespace types

/**
 * Declare type used in Epetra.
 */
using TrilinosScalar = double;


namespace TrilinosWrappers
{
  namespace types
  {
#ifdef DEAL_II_WITH_64BIT_INDICES
    /**
     * Declare type of integer used in the Epetra package of Trilinos.
     */
    using int_type = long long int;
#else
    /**
     * Declare type of integer used in the Epetra package of Trilinos.
     */
    using int_type = int;
#endif
  } // namespace types
} // namespace TrilinosWrappers


// this part of the namespace numbers got moved to the bottom types.h file,
// because otherwise we get a circular inclusion of config.h, types.h, and
// numbers.h
namespace numbers
{
  /**
   * Representation of the largest number that can be put into an unsigned
   * integer. This value is widely used throughout the library as a marker for
   * an invalid unsigned integer value, such as an invalid array index, an
   * invalid array size, and the like.
   */
  static const unsigned int invalid_unsigned_int =
    static_cast<unsigned int>(-1);

  /**
   * Representation of the largest number that can be put into a size_type.
   * This value is used throughout the library as a marker for an invalid
   * size_type value, such as an invalid array index, an invalid array size,
   * and the like. Invalid_size_type is equivalent to invalid_dof_index.
   */
  const types::global_dof_index invalid_size_type =
    static_cast<types::global_dof_index>(-1);

  /**
   * An invalid value for indices of degrees of freedom.
   */
  const types::global_dof_index invalid_dof_index =
    static_cast<types::global_dof_index>(-1);

  /**
   * An invalid value for coarse cell ids. See the glossary
   * entry on @ref GlossCoarseCellId "coarse cell IDs" for more information.
   */
  const types::coarse_cell_id invalid_coarse_cell_id =
    static_cast<types::coarse_cell_id>(-1);

  /**
   * Invalid material_id which we need in several places as a default value.
   * We assume that all material_ids lie in the range [0,
   * invalid_material_id).
   */
  const types::material_id invalid_material_id =
    static_cast<types::material_id>(-1);

  /**
   * Invalid boundary_id which we need in several places as a default value.
   * We assume that all valid boundary_ids lie in the range [0,
   * invalid_boundary_id).
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  const types::boundary_id invalid_boundary_id =
    static_cast<types::boundary_id>(-1);

  /**
   * A boundary indicator number that we reserve for internal faces.  We
   * assume that all valid boundary_ids lie in the range [0,
   * internal_face_boundary_id).
   *
   * This is an indicator that is used internally (by the library) to
   * differentiate between faces that lie at the boundary of the domain and
   * faces that lie in the interior of the domain. You should never try to
   * assign this boundary indicator to anything in user code.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  const types::boundary_id internal_face_boundary_id =
    static_cast<types::boundary_id>(-1);

  /**
   * Invalid manifold_id which we need in several places as a default value.
   * We assume that all valid manifold_ids lie in the range [0,
   * invalid_manifold_id).
   *
   * @deprecated Use flat_manifold_id instead.
   *
   * @see
   * @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  DEAL_II_DEPRECATED
  const types::manifold_id invalid_manifold_id =
    static_cast<types::manifold_id>(-1);

  /**
   * A manifold_id we reserve for the default flat Cartesian manifold.
   *
   * @see
   * @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  const types::manifold_id flat_manifold_id =
    static_cast<types::manifold_id>(-1);

  /**
   * A special id for an invalid subdomain id. This value may not be used as a
   * valid id but is used, for example, for default arguments to indicate a
   * subdomain id that is not to be used.
   *
   * See the
   * @ref GlossSubdomainId "glossary"
   * for more information.
   */
  const types::subdomain_id invalid_subdomain_id =
    static_cast<types::subdomain_id>(-1);

  /**
   * The subdomain id assigned to a cell whose true subdomain id we don't
   * know, for example because it resides on a different processor on a mesh
   * that is kept distributed on many processors. Such cells are called
   * "artificial".
   *
   * See the glossary entries on
   * @ref GlossSubdomainId "subdomain ids"
   * and
   * @ref GlossArtificialCell "artificial cells"
   * as well as the
   * @ref distributed
   * module for more information.
   */
  const types::subdomain_id artificial_subdomain_id =
    static_cast<types::subdomain_id>(-2);
} // namespace numbers

DEAL_II_NAMESPACE_CLOSE

#endif
