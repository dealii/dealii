// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2016 by the deal.II authors
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

#ifndef dealii__cell_id_h
#define dealii__cell_id_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx11/array.h>

#include <vector>
#include <iostream>

#ifdef DEAL_II_WITH_P4EST
#include <deal.II/distributed/p4est_wrappers.h>
#endif

DEAL_II_NAMESPACE_OPEN

template <int, int> class Triangulation;

/**
 * A class to represent a unique ID for a cell in a Triangulation.  This class
 * stores the index of the coarse cell from which a cell is descendant,
 * together with information on how to reach the cell from that coarse cell
 * (i.e., which child index to take on each level when moving from one cell to
 * its children). The important point about this class is that an object of
 * the current class uniquely identifies a cell in triangulation, and it even
 * does so in the context of objects of type
 * parallel::distributed::Triangulation where the local portion of a mesh may
 * not store all cells. For example, the CellId computed for a ghost cell on
 * one processor will be exactly the same as the CellId computed for the very
 * same cell on the processor that actually owns the cell, although the level
 * and index of the iterators pointing to that cell <i>within the
 * triangulation stored on each of the processors</i> may (and in general
 * will) be different. In other words, CellId provides the tool with which it
 * is possible to uniquely identify cells in a parallel triangulation, and
 * consequently makes it possible to exchange data between processors tied to
 * individual cells.
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
   * number of children that can be transferred to 20 in 3D and 30 in 2D
   * (using 2 times 32 bit for storage), a limitation that is identical to
   * the one used by p4est.
   */
  typedef std_cxx11::array<unsigned int, 4> binary_type;

  /**
   * Construct a CellId object with a given @p coarse_cell_id and vector of
   * child indices. @p child_indices is
   * interpreted identical to the member variable with the same name, namely
   * each entry denotes which child to pick from one refinement level to the
   * next, starting with the coarse cell, until we get to the cell represented
   * by the current object. Therefore, each entry should be a number between 0
   * and the number of children of a cell in the current space dimension.
   */
  CellId(const unsigned int coarse_cell_id,
         const std::vector<unsigned char> &child_indices);

  /**
   * Construct a CellId object with a given @p coarse_cell_id and array of
   * child indices provided in @p child_indices. @p child_indices is
   * interpreted identical to the member variable with the same name, namely
   * each entry denotes which child to pick from one refinement level to the
   * next, starting with the coarse cell, until we get to the cell represented
   * by the current object. Therefore, each entry should be a number between 0
   * and the number of children of a cell in the current space dimension.
   * @p child_indices shall have at least @p n_child_indices valid entries.
   */
  CellId(const unsigned int coarse_cell_id,
         const unsigned int n_child_indices,
         const unsigned char *child_indices);

  /**
   * Construct a CellId object with a given binary representation that was
   * previously constructed by CellId::to_binary.
   */
  CellId(const binary_type &binary_representation);

  /**
   * Construct an invalid CellId.
   */
  CellId();

  /**
   * Return a human readable string representation of this CellId.
   */
  std::string to_string() const;

  /**
   * Return a compact and fast binary representation of this CellId.
   */
  template<int dim>
  binary_type to_binary() const;

  /**
   * Return a cell_iterator to the cell represented by this CellId.
   */
  template<int dim, int spacedim>
  typename Triangulation<dim,spacedim>::cell_iterator
  to_cell(const Triangulation<dim,spacedim> &tria) const;

  /**
   * Compare two CellId objects for equality.
   */
  bool operator== (const CellId &other) const;

  /**
   * Compare two CellIds for inequality.
   */
  bool operator!= (const CellId &other) const;

  /**
   * Compare two CellIds with regard to an ordering. The details of this
   * ordering are unspecified except that the operation provides a
   * total ordering among all cells.
   */
  bool operator<(const CellId &other) const;

private:
  /**
   * The number of the coarse cell within whose tree the cell
   * represented by the current object is located.
   */
  unsigned int coarse_cell_id;

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
  std_cxx11::array<char,internal::p4est::functions<2>::max_level> child_indices;
#else
  std_cxx11::array<char,30> child_indices;
#endif

  friend std::istream &operator>> (std::istream &is, CellId &cid);
  friend std::ostream &operator<< (std::ostream &os, const CellId &cid);
};




/**
 * Write a CellId object into a stream.
 */
inline
std::ostream &operator<< (std::ostream &os,
                          const CellId &cid)
{
  os << cid.coarse_cell_id << '_' << cid.n_child_indices << ':';
  for (unsigned int i=0; i<cid.n_child_indices; ++i)
    // write the child indices. because they are between 0 and 2^dim-1, they all
    // just have one digit, so we could write them as integers. it's
    // probably clearer to write them as one-digit characters starting
    // at '0'
    os << static_cast<unsigned char>(cid.child_indices[i] + '0');
  return os;
}



/**
 * Read a CellId object from a stream.
 */
inline
std::istream &operator>> (std::istream &is,
                          CellId &cid)
{
  unsigned int cellid;
  is >> cellid;
  if (is.eof())
    return is;

  cid.coarse_cell_id = cellid;
  char dummy;
  is >> dummy;
  Assert(dummy=='_', ExcMessage("invalid CellId"));
  is >> cid.n_child_indices;
  is >> dummy;
  Assert(dummy==':', ExcMessage("invalid CellId"));

  char value;
  for (unsigned int i=0; i<cid.n_child_indices; ++i)
    {
      // read the one-digit child index (as an integer number) and
      // convert it back into unsigned char
      is >> value;
      cid.child_indices[i] = value-'0';
    }
  return is;
}



inline bool
CellId::operator== (const CellId &other) const
{
  if (this->coarse_cell_id != other.coarse_cell_id)
    return false;
  if (n_child_indices != other.n_child_indices)
    return false;

  for (unsigned int i=0; i<n_child_indices; ++i)
    if (child_indices[i] != other.child_indices[i])
      return false;

  return true;
}



inline bool
CellId::operator!= (const CellId &other) const
{
  return !(*this == other);
}



inline
bool CellId::operator<(const CellId &other) const
{
  if (this->coarse_cell_id != other.coarse_cell_id)
    return this->coarse_cell_id < other.coarse_cell_id;

  unsigned int idx = 0;
  while (idx < n_child_indices)
    {
      if (idx>=other.n_child_indices)
        return false;

      if (child_indices[idx] != other.child_indices[idx])
        return child_indices[idx] < other.child_indices[idx];

      ++idx;
    }

  if (n_child_indices == other.n_child_indices)
    return false;
  return true; // other.id is longer
}

DEAL_II_NAMESPACE_CLOSE

#endif
