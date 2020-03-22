/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

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
 
  ***************************************************************** */

#ifndef MSQ_INTERRUPT_HPP
#define MSQ_INTERRUPT_HPP

#include "MsqError.hpp"

namespace MESQUITE_NS {

/** \brief Class to watch for user-interrupt (SIGINT, ctrl-C)
 *
 * A class to watch for SIGINT.
 * 
 * Creating an instance of the class ensures that the Mesquite
 * handler for SIGINT is registered.  When all instances are
 * destroyed, the Mesquite handler is removed.  The intent is
 * that each interruptable API declare an instance on the stack.
 * This way the handler is automatically unregistered when the 
 * API returns.  For example:
 * <code>
 * void my_api( MsqError& err )
 * {
 *   MsqInterrupt interrupt;
 *   ... //do stuff
 *   return;
 * }
 * </code>
 */
class MsqInterrupt
{
public:

  /**\brief Disable Mesquite's SIGINT handler */
  static void disable( MsqError& err );
  /**\brief Allow Mesquite to register a SIGINT handler */
  static void allow( MsqError& err );
  /**\brief Force Mesquite to register SIGINT handler */
  static void enable( MsqError& err );

  /**\brief Check if an interrupt was seen */
  static bool interrupt()        { return sawInterrupt;  }
  /**\brief Clear the interrupt flag */
  static void clear()            { sawInterrupt = false; }
  /**\brief Set the interrupt flag */
  static void set_interrupt()    { sawInterrupt = true;  } 
  
  /** Constructor, increment instance count. If
    * instance count was zero, register SIGINT handler */
  MsqInterrupt();
  /** Constructor, decrement instance count. If
    * instance count goes to zero, remove SIGINT handler */
  ~MsqInterrupt();
  
  static void set_handler();
  
private:
  
  enum InterruptMode { CATCH, IGNORE, AUTO };
  
  static InterruptMode interruptMode;
  static unsigned instanceCount;
  static bool sawInterrupt;
  
  // Don't allow any of this stuff (make them private)
  void* operator new(size_t size);
  MsqInterrupt( const MsqInterrupt& );
  MsqInterrupt& operator=( const MsqInterrupt& );
};
}

#endif

