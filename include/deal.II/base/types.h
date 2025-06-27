// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_types_h
#define dealii_types_h


#include <deal.II/base/config.h>

#include <cstdint>
#include <type_traits> // make_signed_t


DEAL_II_NAMESPACE_OPEN

/**
 * A namespace in which we declare alias for types used in deal.II, as well
 * as special values for these types.
 */
namespace types
{
  /**
   * The type used to represent face and line orientations.
   *
   * See the
   * @ref GlossCombinedOrientation "glossary"
   * for more information.
   */
  using geometric_orientation = unsigned char;

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
  using global_vertex_index = std::uint64_t;

  /**
   * An identifier that denotes the MPI type associated with
   * types::global_vertex_index.
   *
   * This preprocessor variable is deprecated. Use the variable
   * `Utilities::MPI::mpi_type_id_for_type<types::global_vertex_index>`
   * instead.
   */
#define DEAL_II_VERTEX_INDEX_MPI_TYPE MPI_UINT64_T

  /**
   * The type in which we store the active and future FE indices.
   */
  using fe_index = unsigned short int;

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
  using global_dof_index = std::uint64_t;
#else
  using global_dof_index        = unsigned int;
#endif

  /**
   * Signed version of global_dof_index.
   * This is useful for interacting with Trilinos' Tpetra that only works well
   * with a signed global ordinal type.
   */
#ifdef DEAL_II_WITH_64BIT_INDICES
  using signed_global_dof_index = long long;
#else
  using signed_global_dof_index = int;
#endif

  /**
   * An identifier that denotes the MPI type associated with
   * types::global_dof_index.
   *
   * This preprocessor variable is deprecated. Use the variable
   * `Utilities::MPI::mpi_type_id_for_type<types::global_dof_index>`
   * instead.
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
  using global_cell_index = std::uint64_t;
#else
  using global_cell_index       = unsigned int;
#endif

  /**
   * The type used for coarse-cell ids. See the glossary
   * entry on
   * @ref GlossCoarseCellId "coarse cell IDs"
   * for more information.
   */
  using coarse_cell_id = global_cell_index;

  /**
   * The type used to denote boundary indicators associated with every piece
   * of the boundary of the domain.
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

  /**
   * The type used to denote geometric entity types.
   *
   * @deprecated This type was previously only used in library internals and is
   * deprecated without replacement.
   */
  using geometric_entity_type DEAL_II_DEPRECATED_EARLY = std::uint8_t;
} // namespace types

/**
 * Declare type used in Epetra.
 */
using TrilinosScalar = double;


namespace TrilinosWrappers
{
  namespace types
  {
    /**
     * Declare type of 64 bit integer used in the Epetra package of Trilinos.
     */
    using int64_type = long long int;

#ifdef DEAL_II_WITH_64BIT_INDICES
    /**
     * Declare type of integer used in the Epetra package of Trilinos.
     */
    using int_type = int64_type;
#else
    /**
     * Declare type of integer used in the Epetra package of Trilinos.
     */
    using int_type = int;
#endif
  } // namespace types
} // namespace TrilinosWrappers



namespace numbers
{
  /**
   * Representation of the largest number that can be put into an unsigned
   * integer. This value is widely used throughout the library as a marker for
   * an invalid unsigned integer value, such as an invalid array index, an
   * invalid array size, and the like.
   *
   * This value is an example of an
   * @ref GlossInvalidValue "invalid value".
   * See there for more information.
   */
  constexpr unsigned int invalid_unsigned_int = static_cast<unsigned int>(-1);

  /**
   * Representation of the largest number that can be put into a size_type.
   * This value is used throughout the library as a marker for an invalid
   * size_type value, such as an invalid array index, an invalid array size,
   * and the like. Invalid_size_type is equivalent to invalid_dof_index.
   *
   * This value is an example of an
   * @ref GlossInvalidValue "invalid value".
   * See there for more information.
   */
  constexpr types::global_dof_index invalid_size_type =
    static_cast<types::global_dof_index>(-1);

  /**
   * An invalid value for active and future fe indices.
   *
   * This value is an example of an
   * @ref GlossInvalidValue "invalid value".
   * See there for more information.
   */
  constexpr types::fe_index invalid_fe_index = static_cast<types::fe_index>(-1);

  /**
   * An invalid value for indices of degrees of freedom.
   *
   * This value is an example of an
   * @ref GlossInvalidValue "invalid value".
   * See there for more information.
   */
  constexpr types::global_dof_index invalid_dof_index =
    static_cast<types::global_dof_index>(-1);

  /**
   * An invalid value for coarse cell ids. See the glossary
   * entry on
   * @ref GlossCoarseCellId "coarse cell IDs"
   * for more information.
   *
   * This value is an example of an
   * @ref GlossInvalidValue "invalid value".
   * See there for more information.
   */
  constexpr types::coarse_cell_id invalid_coarse_cell_id =
    static_cast<types::coarse_cell_id>(-1);

  /**
   * Invalid material_id which we need in several places as a default value.
   * We assume that all material_ids lie in the range `[0,
   * invalid_material_id)`.
   *
   * This value is an example of an
   * @ref GlossInvalidValue "invalid value".
   * See there for more information.
   */
  constexpr types::material_id invalid_material_id =
    static_cast<types::material_id>(-1);

  /**
   * Invalid boundary_id which we need in several places as a default value.
   * We assume that all valid boundary_ids lie in the range `[0,
   * invalid_boundary_id)`.
   *
   * This value is an example of an
   * @ref GlossInvalidValue "invalid value".
   * See there for more information.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  constexpr types::boundary_id invalid_boundary_id =
    static_cast<types::boundary_id>(-1);

  /**
   * A boundary indicator number that we reserve for internal faces.  We
   * assume that all valid boundary_ids lie in the range `[0,
   * internal_face_boundary_id)`.
   *
   * This is an indicator that is used internally (by the library) to
   * differentiate between faces that lie at the boundary of the domain and
   * faces that lie in the interior of the domain. You should never try to
   * assign this boundary indicator to anything in user code.
   *
   * This value is an example of an
   * @ref GlossInvalidValue "invalid value".
   * See there for more information.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  constexpr types::boundary_id internal_face_boundary_id =
    static_cast<types::boundary_id>(-1);

  /**
   * A manifold_id we reserve for the default flat Cartesian manifold.
   *
   * This value is an example of an
   * @ref GlossInvalidValue "invalid value".
   * See there for more information.
   *
   * @see
   * @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  constexpr types::manifold_id flat_manifold_id =
    static_cast<types::manifold_id>(-1);

  /**
   * Value indicating that a face or line is in its default orientation.
   *
   * See the
   * @ref GlossCombinedOrientation "glossary"
   * for more information.
   */
  constexpr types::geometric_orientation default_geometric_orientation =
    static_cast<types::geometric_orientation>(0b000);

  /**
   * Value indicating that a line is in the reverse orientation. Since lines can
   * only have two possible orientations, this value and
   * default_geometric_orientation completely encode the possible values for
   * line orientations.
   *
   * See the
   * @ref GlossCombinedOrientation "glossary"
   * for more information.
   */
  constexpr types::geometric_orientation reverse_line_orientation =
    static_cast<types::geometric_orientation>(0b001);

  /**
   * A special id for an invalid subdomain id. This value may not be used as a
   * valid id but is used, for example, for default arguments to indicate a
   * subdomain id that is not to be used.
   *
   * This value is an example of an
   * @ref GlossInvalidValue "invalid value".
   * See there for more information.
   *
   * See the
   * @ref GlossSubdomainId "glossary"
   * for more information.
   */
  constexpr types::subdomain_id invalid_subdomain_id =
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
   * topic for more information.
   *
   * This value is an example of an
   * @ref GlossInvalidValue "invalid value".
   * See there for more information.
   */
  constexpr types::subdomain_id artificial_subdomain_id =
    static_cast<types::subdomain_id>(-2);
} // namespace numbers

DEAL_II_NAMESPACE_CLOSE

#endif
