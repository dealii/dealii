// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_cell_id_h
#define dealii_cell_id_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>

#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

#ifdef DEAL_II_WITH_P4EST
#  include <deal.II/distributed/p4est_wrappers.h>
#endif

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class Triangulation;
#endif

/**
 * A class to represent a unique ID for a cell in a Triangulation. It is
 * returned by `cell->id()` (i.e., CellAccessor::id()) where
 * `cell` is assumed to be a cell iterator.
 *
 * This class stores the index of the coarse cell from which a cell is
 * descendant (or, more specifically, the
 * entry on
 * @ref GlossCoarseCellId "coarse cell IDs"),
 * together with information on how to reach the cell from that coarse cell
 * (i.e., which child index to take on each level of the triangulation when
 * moving from one cell to its children). The important point about this
 * class is that an object of
 * the current class uniquely identifies a cell in triangulation, and it even
 * does so in the context of objects of type
 * parallel::distributed::Triangulation where the local portion of a mesh may
 * not store all cells. For example, the CellId computed for a ghost cell on
 * one processor will be exactly the same as the CellId computed for the very
 * same cell on the processor that actually owns the cell, although the level
 * and index of the iterators pointing to that cell <i>within the
 * triangulation stored on each of the processors</i> may (and in general
 * will) be different. In other words, CellId provides the tool with which it
 * is possible to globally and uniquely identify cells in a parallel
 * triangulation, and consequently makes it possible to exchange, between
 * processors, data tied to individual cells.
 *
 * @note How this data is internally represented is not of importance (and not
 * exposed on purpose).
 */
class CellId
{
public:
  /**
   * A type that is used to encode the CellId data in a compact and fast way
   * (e.g. for MPI transfer to other processes). Note that it limits the
   * number of children that can be transferred to 20 in 3d and 30 in 2d
   * (using 2 times 32 bit for storage), a limitation that is identical to
   * the one used by p4est.
   */
  using binary_type = std::array<unsigned int, 4>;

  /**
   * Construct a CellId object with a given @p coarse_cell_id and vector of
   * child indices. @p child_indices is
   * interpreted identical to the member variable with the same name, namely
   * each entry denotes which child to pick from one refinement level to the
   * next, starting with the coarse cell, until we get to the cell represented
   * by the current object. Therefore, each entry should be a number between 0
   * and the number of children of a cell in the current space dimension (i.e.,
   * GeometryInfo<dim>::max_children_per_cell).
   */
  CellId(const types::coarse_cell_id      coarse_cell_id,
         const std::vector<std::uint8_t> &child_indices);

  /**
   * Construct a CellId object with a given @p coarse_cell_id and array of
   * child indices provided in @p child_indices. @p child_indices is
   * interpreted identical to the member variable with the same name, namely
   * each entry denotes which child to pick from one refinement level to the
   * next, starting with the coarse cell, until we get to the cell represented
   * by the current object. Therefore, each entry should be a number between 0
   * and the number of children of a cell in the current space dimension (i.e.,
   * GeometryInfo<dim>::max_children_per_cell). The array
   * @p child_indices must have at least @p n_child_indices valid entries.
   */
  CellId(const types::coarse_cell_id coarse_cell_id,
         const unsigned int          n_child_indices,
         const std::uint8_t         *child_indices);

  /**
   * Construct a CellId object with a given binary representation that was
   * previously constructed by CellId::to_binary.
   */
  CellId(const binary_type &binary_representation);

  /**
   * Create a CellId from a string with the same format that is produced by
   * to_string().
   */
  explicit CellId(const std::string &string_representation);

  /**
   * Construct an invalid CellId.
   */
  CellId();

  /**
   * Return a human-readable string representation of this CellId.
   *
   * The string returned by this function consists of only ASCII characters
   * and will look, for example, like this: `"0_3:006"`. It *can* be
   * interpreted by humans as saying "This cell originates from the zeroth
   * coarse mesh cell, lives on refinement level 3, and the path from the
   * coarse mesh cell to its children and grand children is given by 006".
   * But it is not *meant* to be interpreted in any meaningful way: It's just
   * a way of representing the internal state of the current object using
   * only ASCII characters in the printable range.
   */
  std::string
  to_string() const;

  /**
   * Return a compact and fast binary representation of this CellId.
   */
  template <int dim>
  binary_type
  to_binary() const;

  /**
   * Compare two CellId objects for equality.
   */
  bool
  operator==(const CellId &other) const;

  /**
   * Compare two CellIds for inequality.
   */
  bool
  operator!=(const CellId &other) const;

  /**
   * Compare two CellIds with regard to an ordering. The details of this
   * ordering are unspecified except that the operation provides a
   * total ordering among all cells.
   */
  bool
  operator<(const CellId &other) const;

  /**
   * Determine if this cell id is the direct parent of the input cell id.
   */
  bool
  is_parent_of(const CellId &other) const;

  /**
   * Determine if this cell id is the ancestor of the input cell id.
   */
  bool
  is_ancestor_of(const CellId &other) const;

  /**
   * Read or write the data of this object to or from a stream for the
   * purpose of serialization using the [BOOST serialization
   * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

  /**
   * Return the id of the coarse cell.
   */
  types::coarse_cell_id
  get_coarse_cell_id() const;

  /**
   * Return a read-only container of integers that denotes which child to pick
   * from one refinement level to the next, starting with the coarse cell, until
   * we get to the cell represented by the current object.
   *
   * The number of elements in this container corresponds to (level-1) of the
   * current cell.
   */
  ArrayView<const std::uint8_t>
  get_child_indices() const;

private:
  /**
   * The number of the coarse cell within whose tree the cell
   * represented by the current object is located.
   */
  types::coarse_cell_id coarse_cell_id;

  /**
   * The number of child indices stored in the child_indices array. This is
   * equivalent to (level-1) of the current cell.
   */
  unsigned int n_child_indices;

  /**
   * An array of integers that denotes which child to pick from one
   * refinement level to the next, starting with the coarse cell,
   * until we get to the cell represented by the current object.
   * Only the first n_child_indices entries are used, but we use a statically
   * allocated array instead of a vector of size n_child_indices to speed up
   * creation of this object. If the given dimensions ever become a limitation
   * the array can be extended.
   */
#ifdef DEAL_II_WITH_P4EST
  std::array<std::uint8_t, internal::p4est::functions<2>::max_level>
    child_indices;
#else
  std::array<std::uint8_t, 30> child_indices;
#endif

  friend std::istream &
  operator>>(std::istream &is, CellId &cid);
  friend std::ostream &
  operator<<(std::ostream &os, const CellId &cid);
};



/**
 * Write a CellId object into a stream.
 */
inline std::ostream &
operator<<(std::ostream &os, const CellId &cid)
{
  os << cid.coarse_cell_id << '_' << cid.n_child_indices << ':';
  for (unsigned int i = 0; i < cid.n_child_indices; ++i)
    // write the child indices. because they are between 0 and 2^dim-1, they all
    // just have one digit, so we could write them as one character
    // objects. it's probably clearer to write them as one-digit characters
    // starting at '0'
    os << static_cast<unsigned char>('0' + cid.child_indices[i]);
  return os;
}



/**
 * Serialization function
 */
template <class Archive>
void
CellId::serialize(Archive &ar, const unsigned int /*version*/)
{
  ar &coarse_cell_id;
  ar &n_child_indices;
  ar &child_indices;
}

/**
 * Read a CellId object from a stream.
 */
inline std::istream &
operator>>(std::istream &is, CellId &cid)
{
  unsigned int cellid;
  is >> cellid;
  if (is.eof())
    return is;

  cid.coarse_cell_id = cellid;
  char dummy;
  is >> dummy;
  Assert(dummy == '_', ExcMessage("invalid CellId"));
  is >> cid.n_child_indices;
  is >> dummy;
  Assert(dummy == ':', ExcMessage("invalid CellId"));

  unsigned char value;
  for (unsigned int i = 0; i < cid.n_child_indices; ++i)
    {
      // read the one-digit child index (as an integer number) and
      // convert it back into unsigned integer type
      is >> value;
      cid.child_indices[i] = value - '0';
    }
  return is;
}



inline bool
CellId::operator==(const CellId &other) const
{
  if (this->coarse_cell_id != other.coarse_cell_id)
    return false;
  if (n_child_indices != other.n_child_indices)
    return false;

  for (unsigned int i = 0; i < n_child_indices; ++i)
    if (child_indices[i] != other.child_indices[i])
      return false;

  return true;
}



inline bool
CellId::operator!=(const CellId &other) const
{
  return !(*this == other);
}



inline bool
CellId::operator<(const CellId &other) const
{
  if (this->coarse_cell_id != other.coarse_cell_id)
    return this->coarse_cell_id < other.coarse_cell_id;

  unsigned int idx = 0;
  while (idx < n_child_indices)
    {
      if (idx >= other.n_child_indices)
        return false;

      if (child_indices[idx] != other.child_indices[idx])
        return child_indices[idx] < other.child_indices[idx];

      ++idx;
    }

  if (n_child_indices == other.n_child_indices)
    return false;
  return true; // other.id is longer
}



inline bool
CellId::is_parent_of(const CellId &other) const
{
  if (this->coarse_cell_id != other.coarse_cell_id)
    return false;

  if (n_child_indices + 1 != other.n_child_indices)
    return false;

  for (unsigned int idx = 0; idx < n_child_indices; ++idx)
    if (child_indices[idx] != other.child_indices[idx])
      return false;

  return true; // other.id is longer
}



inline bool
CellId::is_ancestor_of(const CellId &other) const
{
  if (this->coarse_cell_id != other.coarse_cell_id)
    return false;

  if (n_child_indices >= other.n_child_indices)
    return false;

  for (unsigned int idx = 0; idx < n_child_indices; ++idx)
    if (child_indices[idx] != other.child_indices[idx])
      return false;

  return true; // other.id is longer
}



inline types::coarse_cell_id
CellId::get_coarse_cell_id() const
{
  return coarse_cell_id;
}



inline ArrayView<const std::uint8_t>
CellId::get_child_indices() const
{
  return {child_indices.data(), n_child_indices};
}


DEAL_II_NAMESPACE_CLOSE

#endif
