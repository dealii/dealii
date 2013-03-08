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

#ifndef __TBB_machine_H
#error Do not include this file directly; include tbb_machine.h instead
#endif

#include <sched.h>
#define __TBB_Yield()  sched_yield()


// __TBB_DetectNumberOfWorkers

#include <sys/types.h>
#include <sys/sysctl.h>

static inline int __TBB_macos_available_cpu() {
    int name[2] = {CTL_HW, HW_AVAILCPU};
    int ncpu;
    size_t size = sizeof(ncpu);
    sysctl( name, 2, &ncpu, &size, NULL, 0 );
    return ncpu;
}

#define __TBB_DetectNumberOfWorkers() __TBB_macos_available_cpu()


#ifndef __TBB_WORDSIZE
#define __TBB_WORDSIZE 4
#endif

#ifndef __TBB_BIG_ENDIAN
#if __BIG_ENDIAN__
#define __TBB_BIG_ENDIAN 1
#else
#define __TBB_BIG_ENDIAN 0
#endif
#endif


#if !defined(__TBB_CompareAndSwap4) || !defined(__TBB_CompareAndSwap8)

// Implementation of atomic operations based on OS provided primitives
#include <libkern/OSAtomic.h>

#define __TBB_release_consistency_helper() OSMemoryBarrier()
#define __TBB_full_memory_fence()          OSMemoryBarrier()

static inline int32_t __TBB_macos_cmpswp4(volatile void *ptr, int32_t value, int32_t comparand)
{
    __TBB_ASSERT( !((uintptr_t)ptr&0x3), "address not properly aligned for Mac OS atomics");
    int32_t* address = (int32_t*)ptr;
    while( !OSAtomicCompareAndSwap32Barrier(comparand, value, address) ){
        int32_t snapshot = *address;
        if( snapshot!=comparand ) return snapshot;
    }
    return comparand;
}

static inline int64_t __TBB_macos_cmpswp8(volatile void *ptr, int64_t value, int64_t comparand)
{
    __TBB_ASSERT( !((uintptr_t)ptr&0x7), "address not properly aligned for Mac OS atomics");
    int64_t* address = (int64_t*)ptr;
    while( !OSAtomicCompareAndSwap64Barrier(comparand, value, address) ){
#if __TBB_WORDSIZE==8
        int64_t snapshot = *address;
#else
        int64_t snapshot = OSAtomicAdd64( 0, address );
#endif
        if( snapshot!=comparand ) return snapshot;
    }
    return comparand;
}

#define __TBB_CompareAndSwap4(P,V,C) __TBB_macos_cmpswp4(P,V,C)
#define __TBB_CompareAndSwap8(P,V,C) __TBB_macos_cmpswp8(P,V,C)

static inline int32_t __TBB_macos_fetchadd4(volatile void *ptr, int32_t addend)
{
    __TBB_ASSERT( !((uintptr_t)ptr&0x3), "address not properly aligned for Mac OS atomics");
    return OSAtomicAdd32Barrier(addend, (int32_t*)ptr) - addend;
}

static inline int64_t __TBB_macos_fetchadd8(volatile void *ptr, int64_t addend)
{
    __TBB_ASSERT( !((uintptr_t)ptr&0x7), "address not properly aligned for Mac OS atomics");
    return OSAtomicAdd64Barrier(addend, (int64_t*)ptr) - addend;
}

#define __TBB_FetchAndAdd4(P,V) __TBB_macos_fetchadd4(P,V)
#define __TBB_FetchAndAdd8(P,V) __TBB_macos_fetchadd8(P,V)

#if __TBB_WORDSIZE==4
#define __TBB_CompareAndSwapW(P,V,C) __TBB_CompareAndSwap4(P,V,C)
#define __TBB_FetchAndAddW(P,V) __TBB_FetchAndAdd4(P,V)
#else
#define __TBB_CompareAndSwapW(P,V,C) __TBB_CompareAndSwap8(P,V,C)
#define __TBB_FetchAndAddW(P,V) __TBB_FetchAndAdd8(P,V)
#endif

#endif /* !defined(__TBB_CompareAndSwap4) || !defined(__TBB_CompareAndSwap8) */
