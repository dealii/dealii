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


/** \file Bits.hpp
 *  \brief Functions implementing misc. bit-wise algorithms
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_BITS_HPP
#define MSQ_BITS_HPP

#include "Mesquite.hpp"

namespace MESQUITE_NS {

/**\brief Count number of 1-bits in an unsigned integer
 *
 * This operation is typically referred to as 'popcount' (short for
 * "population count".)  For example the GCC builtin "__builtin_popcount".
 * It is also sometimes referred to as "sideways addition" (e.g. Knuth
 * volume 4.)
 */
inline int popcount( unsigned int x ) {
#if 0  // April '08 bug report says gcc builtin not so good unless Intel SSE4
//#ifdef __GNUC__
  return __builtin_popcount( x );
#else
    // Use "parallel" algorithm
//  if (sizeof(x) == 4) {
    x = (x & 0x55555555u) + ((x >> 1) & 0x55555555u);
    x = (x & 0x33333333u) + ((x >> 2) & 0x33333333u);
    x = (x & 0x0f0f0f0fu) + ((x >> 4) & 0x0f0f0f0fu);
    return (x * 0x1010101u) >> 24; // x % 255
//  }
//  else {
//    x = (x & 0x5555555555555555ULL) + ((x >> 1) & 0x5555555555555555ULL);
//    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
//    x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
//    return (x * 0x0101010101010101ULL) >> 56;
//  }    
#endif
}

} // namespace Mesquite

#endif
