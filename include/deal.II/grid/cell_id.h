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

#include <vector>
#include <iostream>


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
 *
 * @todo Does it make sense to implement a more efficient representation
 * (internally and/or as a string)? If yes, something like a 64bit int as in
 * p4est would be a good option.
 */
class CellId
{
public:
  /**
   * Construct a CellId object with a given @p coarse_cell_index and list of child indices.
   */
  CellId(const unsigned int coarse_cell_id,
         const std::vector<unsigned char> &child_indices);

  /**
   * Construct an invalid CellId.
   */
  CellId();

  /**
   * Return a string representation of this CellId.
   */
  std::string to_string() const;

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
   * A list of integers that denote which child to pick from one
   * refinement level to the next, starting with the coarse cell,
   * until we get to the cell represented by the current object.
   */
  std::vector<unsigned char> child_indices;

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
  os << cid.coarse_cell_id << '_' << cid.child_indices.size() << ':';
  for (unsigned int i=0; i<cid.child_indices.size(); ++i)
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
  unsigned int idsize;
  is >> idsize;
  is >> dummy;
  Assert(dummy==':', ExcMessage("invalid CellId"));

  char value;
  cid.child_indices.clear();
  for (unsigned int i=0; i<idsize; ++i)
    {
      // read the one-digit child index (as an integer number) and
      // convert it back into unsigned char
      is >> value;
      cid.child_indices.push_back(value-'0');
    }
  return is;
}


inline
CellId::CellId(const unsigned int coarse_cell_id,
               const std::vector<unsigned char> &id)
  :
  coarse_cell_id(coarse_cell_id),
  child_indices(id)
{}


inline
CellId::CellId()
  :
  coarse_cell_id(static_cast<unsigned int>(-1))
{}



inline bool
CellId::operator== (const CellId &other) const
{
  if (this->coarse_cell_id != other.coarse_cell_id)
    return false;
  return child_indices == other.child_indices;
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
  while (idx < child_indices.size())
    {
      if (idx>=other.child_indices.size())
        return false;

      if (child_indices[idx] != other.child_indices[idx])
        return child_indices[idx] < other.child_indices[idx];

      ++idx;
    }

  if (child_indices.size() == other.child_indices.size())
    return false;
  return true; // other.id is longer
}

DEAL_II_NAMESPACE_CLOSE

#endif
