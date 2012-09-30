/*
    Copyright 2005-2010 Intel Corporation.  All Rights Reserved.

    This file is part of Threading Building Blocks.

    Threading Building Blocks is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    version 2 as published by the Free Software Foundation.

    Threading Building Blocks is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Threading Building Blocks; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    As a special exception, you may use this file as part of a free software
    library without restriction.  Specifically, if other files instantiate
    templates or use macros or inline functions from this file, or you compile
    this file and link it with other files to produce an executable, this
    file does not by itself cause the resulting executable to be covered by
    the GNU General Public License.  This exception does not however
    invalidate any other reasons why the executable file might be covered by
    the GNU General Public License.
*/

#ifndef _TBB_tbb_misc_H
#define _TBB_tbb_misc_H

#include "tbb/tbb_stddef.h"
#include "tbb/tbb_machine.h"

#if _WIN32||_WIN64
#include "tbb/machine/windows_api.h"
#elif __linux__
#include <sys/sysinfo.h>
#define __TBB_DetectNumberOfWorkers() get_nprocs()
#elif defined(__sun)
#include <sys/sysinfo.h>
#include <unistd.h>
#elif defined(__NetBSD__) || defined(__FreeBSD__) || defined(_AIX)
#include <unistd.h>
#endif

namespace tbb {
namespace internal {

const size_t MByte = 1<<20;

#if !defined(__TBB_WORDSIZE)
    const size_t ThreadStackSize = 1*MByte;
#elif __TBB_WORDSIZE<=4
    const size_t ThreadStackSize = 2*MByte;
#else
    const size_t ThreadStackSize = 4*MByte;
#endif

#if defined(__TBB_DetectNumberOfWorkers) // covers Linux, Mac OS*, and other platforms

static inline int DetectNumberOfWorkers() {
    int n = __TBB_DetectNumberOfWorkers(); 
    return n>0? n: 1; // Fail safety strap
}

#else /* !__TBB_DetectNumberOfWorkers */

#if _WIN32||_WIN64

static inline int DetectNumberOfWorkers() {
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    return static_cast<int>(si.dwNumberOfProcessors);
}

#elif defined(_SC_NPROCESSORS_ONLN)

static inline int DetectNumberOfWorkers() {
    int number_of_workers = sysconf(_SC_NPROCESSORS_ONLN);
    return number_of_workers>0? number_of_workers: 1;
}

#else
#error DetectNumberOfWorkers: Method to detect the number of available CPUs is unknown
#endif /* os kind */

#endif /* !__TBB_DetectNumberOfWorkers */

//! Throws std::runtime_error with what() returning error_code description prefixed with aux_info
void handle_win_error( int error_code );

//! True if environment variable with given name is set and not 0; otherwise false.
bool GetBoolEnvironmentVariable( const char * name );

//! Print TBB version information on stderr
void PrintVersion();

//! Print extra TBB version information on stderr
void PrintExtraVersionInfo( const char* category, const char* description );

//! A callback routine to print RML version information on stderr
void PrintRMLVersionInfo( void* arg, const char* server_info );

// For TBB compilation only; not to be used in public headers
#if defined(min) || defined(max)
#undef min
#undef max
#endif

//! Utility template function returning lesser of the two values.
/** Provided here to avoid including not strict safe <algorithm>.\n
    In case operands cause signed/unsigned or size mismatch warnings it is caller's
    responsibility to do the appropriate cast before calling the function. **/
template<typename T1, typename T2>
T1 min ( const T1& val1, const T2& val2 ) {
    return val1 < val2 ? val1 : val2;
}

//! Utility template function returning greater of the two values.
/** Provided here to avoid including not strict safe <algorithm>.\n
    In case operands cause signed/unsigned or size mismatch warnings it is caller's
    responsibility to do the appropriate cast before calling the function. **/
template<typename T1, typename T2>
T1 max ( const T1& val1, const T2& val2 ) {
    return val1 < val2 ? val2 : val1;
}

//------------------------------------------------------------------------
// FastRandom
//------------------------------------------------------------------------

/** Defined in tbb_main.cpp **/
unsigned GetPrime ( unsigned seed );

//! A fast random number generator.
/** Uses linear congruential method. */
class FastRandom {
    unsigned x, a;
public:
    //! Get a random number.
    unsigned short get() {
        return get(x);
    }
    //! Get a random number for the given seed; update the seed for next use.
    unsigned short get( unsigned& seed ) {
        unsigned short r = (unsigned short)(seed>>16);
        seed = seed*a+1;
        return r;
    }
    //! Construct a random number generator.
    FastRandom( unsigned seed ) {
        x = seed;
        a = GetPrime( seed );
    }
};

} // namespace internal
} // namespace tbb

#endif /* _TBB_tbb_misc_H */
