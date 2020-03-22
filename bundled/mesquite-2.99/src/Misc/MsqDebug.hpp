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

#ifndef MSQ_DEBUG_HPP
#define MSQ_DEBUG_HPP


#include "Mesquite.hpp"

#include <vector>
#include <iostream>
#include <cstdio>

namespace MESQUITE_NS {

/**
 *\defgroup debug Mesquite debug output
 */
/*@{*/

/** \def MSQ_ENABLE_DEBUG
 *\brief Enable debug output and optionally set initial debug flags
 *
 *If this preprocessor constant is undefined, all debug output will
 *be disabled.  It is expected that this will not be defined for
 *release builds of Mesquite.
 *
 *If MSQ_ENABLE_DEBUG is defined, it may optionally be set to a
 *comma-separated list of positive, non-zero integers.  This list
 *will be interpreted as a list of debug flags that should be 
 *activated automatically.  
 *
 *To enable debug output but leave all debug flags deactivated:
 *#define  MSQ_ENABLE_DEBUG
 *
 *To enable debug output and activated the first three debug flags:
 *#define  MSQ_ENABLE_DEBUG 1,2,3
 */

/**
 *\class MsqDebug
 *\brief Run-time activation/deactivation of debug flags
 *\author Jason Kraftcheck
 *\date 2004-10-18
 *
 *Wrap static functions for managing debug flags and associated
 *output streams.  Output is expected do be done using the 
 *provided macros MSQ_DBGOUT(), MSQ_PRINT(), and MSQ_DBG().
 *
 *The default output stream for all flags is cout.
 */
class MsqDebug
{
  public:
    
      // some pre-defined meanings of debug flags
    enum {
      WARN = 1,
      INFO = 2
    };
  
      /**\brief Enable a debug flag */
    static void enable( unsigned flag )          { set( flag, true );  }
      /**\brief Disable a debug flag */
    static void disable( unsigned flag )         { set( flag, false ); }
      /**\brief Set a debug flag */
    static void set( unsigned flag, bool state );
      /**\brief Check a debug flag */
    static bool get( unsigned flag );
    
      /**\brief Disable all debug streams */
    static void disable_all();
  
      /**\brief Get the output stream to be used for a given debug flag */
    static std::ostream& get_stream( unsigned flag );
      /**\brief Set the output stream to be used for a given debug flag */
    static void set_stream( unsigned flag, std::ostream& stream );
    
      // Work around limitations of preprocessor macros.
      // You probably don't want to use this directly.  See
      // MSQ_PRINT macro below.
    class FormatPrinter {
      public:
        FormatPrinter( unsigned flag ) : myFlag(flag) {}
        void print( const char* format, ... ) const
          #ifdef __GNUC__
           __attribute__ ((format (printf, 2, 3)))
          #endif
          ;
        const unsigned myFlag;
    };
    
      // Static initialize function (declare a static instance of this
      // such that the constructor can be used as a function to initialize
      // static data.)
    class InitializeFlags {
      public:
        InitializeFlags();
    };
      
  private:
  
    static std::vector<std::ostream*> streams;
    static std::vector<bool> flags;
    static InitializeFlags init;
};

/** \brief Check if a debug flag is activated - evaluates to a bool.
 */
#ifdef MSQ_ENABLE_DEBUG
#  define MSQ_DBG(flag) Mesquite::MsqDebug::get(flag)
#else
#  define MSQ_DBG(flag) false
#endif

/**
 *\brief Check debug flag and return ostream associated with flag.
 *
 * Evaluates to a conditional calling of an ostream depending on the 
 * debug flag.  Example:
 *   MSQ_DBGOUT(f) << "Debug flag " << f << " is activated" << endl;
 */
#define MSQ_DBGOUT(flag) if (MSQ_DBG(flag)) Mesquite::MsqDebug::get_stream(flag)

/**
 *\brief Check debug flag and print printf-style formatted output.
 *
 * Evaluatues to a conditional print depending on the debug flag.
 * Example:
 *  MSQ_PRINT(f)("Debug flag %d is activated", f);
 */
#define MSQ_PRINT(flag)  if (MSQ_DBG(flag)) Mesquite::MsqDebug::FormatPrinter(flag).print

/*@}*/

} // namespace Mesquite

#endif
