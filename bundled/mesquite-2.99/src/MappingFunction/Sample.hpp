/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file Sample.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_SAMPLE_HPP
#define MSQ_SAMPLE_HPP

#include "Mesquite.hpp"
#include <iosfwd>
#include <stddef.h> // for size_t

namespace MESQUITE_NS {

//! Define a logical location within an element at which the
//! element will be "sampled" for the purpose of evaluating a
//! quality metric.  For example, specify a location at which
//! to calculate the Jacobian of the mapping function and
//! use it to evaluate a TMP quality metric.
//!
//! Logical sample points are currently limited to corners, mid-sides
//! (region elements have both face and edge "sides") and mid-element.
struct MESQUITE_EXPORT Sample {
  //!\brief The "dimension" of the sub-entity.  
  //!
  //! Possible values are [0-3]:
  //! - 0 : sample at the location of a corner vertex
  //! - 1 : sample at the center of an edge
  //! - 2 : sample at the center of a surface element or the center
  //!       of one face of a volume element
  //! - 3 : sample at the center of a volume element
  unsigned short dimension;
  
  //!\brief Canonical number in ITAPS ordering for element corners/sides.
  //!
  //! The dimension specifies which set of sub-entities that the sample
  //! my be on.  The number specifies which entity within that subset.
  //! For example, if dimension == 1 then the sample point is at the
  //! center of one of the element edges.  The number then specifies
  //! which edge of the element, where edges are ordered according to
  //! the ITAPS ordering.
  //!
  //! For mid-region or mid-face of a surface element, the number
  //! should be zero.
  unsigned short number;
  
  //! Constants used for a packed (minimum number of bits) representation
  //! of a sample.
  enum {
    /** Number of bits used to store the dimension of an element 'side' */
    SIDE_DIMENSION_BITS = 2,
    /** Number of bits used to store the index of an element 'side' of a specific dimension */
    SIDE_NUMBER_BITS = 4,
    NUMBER_PACKED_BITS = SIDE_DIMENSION_BITS + SIDE_NUMBER_BITS,
    /** Number of distinct side dimension values that will fit
     *  in a sample value (one greater than the largest dimension) */
    NUM_SAMPLE_SIDE_DIM = 1u << SIDE_DIMENSION_BITS,
    /** Number of distinct side index values that will fit
     *  in a sample value (one greater than the largest side number) */
    NUM_SAMPLE_SIDE_NUM = 1u << SIDE_NUMBER_BITS,
    /** Mask to remove side dimension bits from sample number */
    SIDE_NUMBER_MASK = NUM_SAMPLE_SIDE_NUM - 1,
    /** Mask to remove side dimension bits from sample number */
    SIDE_DIMENSON_MASK = NUM_SAMPLE_SIDE_DIM - 1
  };
  
  //! Return packed representation of this sample.
  inline size_t pack() const
    { return (((size_t)dimension) << SIDE_NUMBER_BITS) | number; }
  //! Set this sample to the values encoded in a packed sample representation.
  inline void unpack( size_t packed )
    { dimension = (packed >> SIDE_NUMBER_BITS) & SIDE_DIMENSON_MASK;
      number = packed & SIDE_NUMBER_MASK;
    }
  //! Initialization constructor
  Sample( unsigned dim, unsigned num ) : dimension(dim), number(num) {}
  //! Initialize from packed representation
  explicit Sample( size_t packed ) { unpack(packed); }
  //! Default constructor must be explicitly included if we provide other constructors
  //! Do nothing (don't waste time initilazing to zero or something, I'd rather
  //! be able to catch the use of uninitialized values using a memory checker
  //! anyway if I make such a mistake.)
  Sample() {}
  
  bool operator==( const Sample& other ) const 
    { return pack() == other.pack(); }
  bool operator!=( const Sample& other ) const 
    { return pack() != other.pack(); }
  bool operator<( const Sample& other ) const 
    { return pack() < other.pack(); }
  bool operator>( const Sample& other ) const 
    { return pack() > other.pack(); }
  bool operator<=( const Sample& other ) const 
    { return pack() <= other.pack(); }
  bool operator>=( const Sample& other ) const 
    { return pack() >= other.pack(); }
};

std::ostream& operator <<( std::ostream& str, Sample s );

} // namespace Mesquite

#endif
