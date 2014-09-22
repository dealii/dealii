// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2014 by the deal.II authors
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

#ifndef __deal2__cell_id_h
#define __deal2__cell_id_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN


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
 * (internally and/or as a string)? If yes, something like a 64bit int
 * as in p4est would be a good option.
 */
class CellId
{
public:
  /**
   * construct CellId with a given coarse_cell_index and list of child indices
   */
  explicit CellId(unsigned int coarse_cell_id_, std::vector<unsigned char> id_)
    : coarse_cell_id(coarse_cell_id_), id(id_)
  {}

  /**
   * construct an empty CellId.
   */
  CellId()
    : coarse_cell_id(-1)
  {}

  /**
   * compare two CellIds
   */
  bool operator== (const CellId &other) const;

  /**
   * compare two CellIds
   */
  bool operator!= (const CellId &other) const;

  /**
   * compare two CellIds
   */
  bool operator<(const CellId &other) const;

  friend std::istream &operator>> (std::istream &is, CellId &cid);
  friend std::ostream &operator<< (std::ostream &os, const CellId &cid);
private:
  unsigned int coarse_cell_id;
  std::vector<unsigned char> id;
};

/**
 * output CellId into a stream
 */
inline std::ostream &operator<< (std::ostream &os, const CellId &cid)
{
  os << cid.coarse_cell_id << '_' << cid.id.size() << ':';
  for (unsigned int i=0; i<cid.id.size(); ++i)
    os << static_cast<int>(cid.id[i]);
  return os;
}

/**
 * read CellId from a stream
 */
inline std::istream &operator>> (std::istream &is, CellId &cid)
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
  cid.id.clear();
  for (unsigned int i=0; i<idsize; ++i)
    {
      is >> value;
      cid.id.push_back(value-'0');
    }
  return is;
}

inline bool
CellId::operator== (const CellId &other) const
{
  if (this->coarse_cell_id != other.coarse_cell_id)
    return false;
  return id == other.id;
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
  while (idx < id.size())
    {
      if (idx>=other.id.size())
        return false;

      if (id[idx] != other.id[idx])
        return id[idx] < other.id[idx];

      ++idx;
    }

  if (id.size() == other.id.size())
    return false;
  return true; // other.id is longer
}

DEAL_II_NAMESPACE_CLOSE

#endif
