// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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
 * A class to represent a unique ID for a cell in a Triangulation.
 * This class stores the index of the coarse cell together with the
 * information on how to reach the cell from that coarse cell (which
 * child index to take on each level). The internal representation
 * is not exposed on purpose.
 *
 * TODO: does it make sense to implement a more efficient representation
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

/**
 *
 */
inline bool
CellId::operator!= (const CellId &other) const
{
  return !(*this == other);
}

DEAL_II_NAMESPACE_CLOSE

#endif
