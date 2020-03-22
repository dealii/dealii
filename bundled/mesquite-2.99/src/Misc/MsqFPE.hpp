/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2005 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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
 
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */

#ifndef MSQ_FPE_HPP
#define MSQ_FPE_HPP

#include <sys/types.h>
#include "Mesquite.hpp"

namespace MESQUITE_NS {

/**\brief Utility class used by InstructionQueue SIGFPE option
 *
 * This is a simple utility class for enabling floating point
 * exceptions.  It provides two functionalities.  The first,
 * implemented in the static methods, is a platform-independent
 * mechanism for modifying the the FPE state.  The second,
 * implemented in the constructor/destructor, is a utlity for
 * setting and resetting the FPE state.  The FPE state is set
 * when the object is created and reset when the object is 
 * destroted.  The intention is that an instance of this object
 * be declared on the stack such that when the instantiating function
 * returns, the destructor is automatically invoked, restoring the
 * original state.
 */
class MsqFPE
{
public:

    /**\brief Set FPE state
     *
     * If <code>enabled == true</code>, floating point exceptions
     * are enabled by the constructor and reset by the destructor.
     * If <code>enabled == false</code>, nothing is done.
     */
  MsqFPE( bool enabled );
  
    /**\brief Restore FPE state */
  ~MsqFPE();
  
    /**\brief Check if FPE state manipulation is supported on this platform */
  static bool fpe_trap_supported();
    /**\return An integer representing the current FPE flags */
  static int get_current_fpe_state();
    /**\return Set the FPE flags on the processor */
  static void set_current_fpe_state(int state);
    /**\return Enable trapping of INVALID, DIVBYZERO, and OVERFLOW */
  static void enable_trap_fpe();
  
private:

    /**\brief dummy declaration preventing heap allocation */
  void* operator new(size_t);
    /**\brief dummy declaration preventing default copy constructor */
  MsqFPE( const MsqFPE& );
    /**\brief dummy declaration preventing default assignment operator */
  MsqFPE& operator=( const MsqFPE& );

    /** Saved constructor argument for use in destructor */
  bool isEnabled;
    /** Saved FPE state for use in destructor */
  int prevState;
};

} // namespace Mesquite

#endif

