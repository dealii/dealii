/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file SimpleStats.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_SIMPLE_STATS_HPP
#define MSQ_SIMPLE_STATS_HPP

#include "Mesquite.hpp"

namespace MESQUITE_NS {

/**\brief Accumulate various statistics for a list of discrete values */
class SimpleStats 
{
private:
  double valueSum, valueSqrSum, valueMin, valueMax;
  unsigned long valueCount;

public:
  SimpleStats();
  
  /**\brief minimum value */
  double minimum()  const { return valueMin; }
  /**\brief maximum value */
  double maximum()  const { return valueMax; }
  /**\brief algebraic mean of values */
  double average()  const { return valueSum/valueCount; }
  /**\brief root mean squared of values */
  double rms()      const { return sqrt(valueSqrSum/valueCount); }
  /**\brief variance of values */
  double variance() const { return valueSqrSum/valueCount - average()*average(); }
  /**\brief standard deviation of values */
  double standard_deviation() const { return sqrt(fabs(variance())); }

  /**\brief incorporate another value into statistics */
  void add_value( double value ) {
    valueSum += value;
    valueSqrSum += value*value;
    if (value < valueMin)
      valueMin = value;
    if (value > valueMax)
      valueMax = value;
    ++valueCount;
  }

  /**\brief incorporate another value into statistics */
  void add_squared( double value_squared ) {
    double value = sqrt(value_squared);
    valueSum += value;
    valueSqrSum += value_squared;
    if (value < valueMin)
      valueMin = value;
    if (value > valueMax)
      valueMax = value;
    ++valueCount;
  }
  
  void clear();
  
  bool empty() { return 0ul == valueCount; }
};


} // namespace MESQUITE_NS

#endif
