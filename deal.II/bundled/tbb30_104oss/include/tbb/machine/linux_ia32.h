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

#include <stdint.h>
#include <unistd.h>

#define __TBB_WORDSIZE 4
#define __TBB_BIG_ENDIAN 0

#define __TBB_release_consistency_helper() __asm__ __volatile__("": : :"memory")
#define __TBB_full_memory_fence() __asm__ __volatile__("mfence": : :"memory")

#if __TBB_ICC_ASM_VOLATILE_BROKEN
#define __TBB_VOLATILE
#else
#define __TBB_VOLATILE volatile
#endif

#define __MACHINE_DECL_ATOMICS(S,T,X) \
static inline T __TBB_machine_cmpswp##S (volatile void *ptr, T value, T comparand )  \
{                                                                                    \
    T result;                                                                        \
                                                                                     \
    __asm__ __volatile__("lock\ncmpxchg" X " %2,%1"                                  \
                          : "=a"(result), "=m"(*(__TBB_VOLATILE T*)ptr)              \
                          : "q"(value), "0"(comparand), "m"(*(__TBB_VOLATILE T*)ptr) \
                          : "memory");                                               \
    return result;                                                                   \
}                                                                                    \
                                                                                     \
static inline T __TBB_machine_fetchadd##S(volatile void *ptr, T addend)              \
{                                                                                    \
    T result;                                                                        \
    __asm__ __volatile__("lock\nxadd" X " %0,%1"                                     \
                          : "=r"(result), "=m"(*(__TBB_VOLATILE T*)ptr)              \
                          : "0"(addend), "m"(*(__TBB_VOLATILE T*)ptr)                \
                          : "memory");                                               \
    return result;                                                                   \
}                                                                                    \
                                                                                     \
static inline  T __TBB_machine_fetchstore##S(volatile void *ptr, T value)            \
{                                                                                    \
    T result;                                                                        \
    __asm__ __volatile__("lock\nxchg" X " %0,%1"                                     \
                          : "=r"(result), "=m"(*(__TBB_VOLATILE T*)ptr)              \
                          : "0"(value), "m"(*(__TBB_VOLATILE T*)ptr)                 \
                          : "memory");                                               \
    return result;                                                                   \
}                                                                                    \
                                                                                     
__MACHINE_DECL_ATOMICS(1,int8_t,"")
__MACHINE_DECL_ATOMICS(2,int16_t,"")
__MACHINE_DECL_ATOMICS(4,int32_t,"l")

static inline int64_t __TBB_machine_cmpswp8 (volatile void *ptr, int64_t value, int64_t comparand )
{
    int64_t result;
#if __PIC__ 
    /* compiling position-independent code */
    // EBX register preserved for compliance with position-independent code rules on IA32
    __asm__ __volatile__ (
            "pushl %%ebx\n\t"
            "movl  (%%ecx),%%ebx\n\t"
            "movl  4(%%ecx),%%ecx\n\t"
            "lock\n\t cmpxchg8b %1\n\t"
            "popl  %%ebx"
             : "=A"(result), "=m"(*(int64_t *)ptr)
             : "m"(*(int64_t *)ptr)
             , "0"(comparand)
             , "c"(&value)
             : "memory", "esp"
#if __INTEL_COMPILER
             ,"ebx"
#endif
    );
#else /* !__PIC__ */
    union {
        int64_t i64;
        int32_t i32[2];
    };
    i64 = value;
    __asm__ __volatile__ (
            "lock\n\t cmpxchg8b %1\n\t"
             : "=A"(result), "=m"(*(__TBB_VOLATILE int64_t *)ptr)
             : "m"(*(__TBB_VOLATILE int64_t *)ptr)
             , "0"(comparand)
             , "b"(i32[0]), "c"(i32[1])
             : "memory"
    );
#endif /* __PIC__ */
    return result;
}

static inline int32_t __TBB_machine_lg( uint32_t x ) {
    int32_t j;
    __asm__ ("bsr %1,%0" : "=r"(j) : "r"(x));
    return j;
}

static inline void __TBB_machine_or( volatile void *ptr, uint32_t addend ) {
    __asm__ __volatile__("lock\norl %1,%0" : "=m"(*(__TBB_VOLATILE uint32_t *)ptr) : "r"(addend), "m"(*(__TBB_VOLATILE uint32_t *)ptr) : "memory");
}

static inline void __TBB_machine_and( volatile void *ptr, uint32_t addend ) {
    __asm__ __volatile__("lock\nandl %1,%0" : "=m"(*(__TBB_VOLATILE uint32_t *)ptr) : "r"(addend), "m"(*(__TBB_VOLATILE uint32_t *)ptr) : "memory");
}

static inline void __TBB_machine_pause( int32_t delay ) {
    for (int32_t i = 0; i < delay; i++) {
       __asm__ __volatile__("pause;");
    }
    return;
}   

static inline int64_t __TBB_machine_load8 (const volatile void *ptr) {
    int64_t result;
    if( ((uint32_t)ptr&7u)==0 ) {
        // Aligned load
        __asm__ __volatile__ ( "fildq %1\n\t"
                               "fistpq %0" :  "=m"(result) : "m"(*(const __TBB_VOLATILE uint64_t*)ptr) : "memory" );
    } else {
        // Unaligned load
        result = __TBB_machine_cmpswp8(const_cast<void*>(ptr),0,0);
    }
    return result;
}

//! Handles misaligned 8-byte store
/** Defined in tbb_misc.cpp */
extern "C" void __TBB_machine_store8_slow( volatile void *ptr, int64_t value );
extern "C" void __TBB_machine_store8_slow_perf_warning( volatile void *ptr );

static inline void __TBB_machine_store8(volatile void *ptr, int64_t value) {
    if( ((uint32_t)ptr&7u)==0 ) {
        // Aligned store
        __asm__ __volatile__ ( "fildq %1\n\t"
                               "fistpq %0" :  "=m"(*(__TBB_VOLATILE int64_t*)ptr) : "m"(value) : "memory" );
    } else {
        // Unaligned store
#if TBB_USE_PERFORMANCE_WARNINGS
        __TBB_machine_store8_slow_perf_warning(ptr);
#endif /* TBB_USE_PERFORMANCE_WARNINGS */
        __TBB_machine_store8_slow(ptr,value);
    }
}
 
// Machine specific atomic operations

#define __TBB_CompareAndSwap1(P,V,C) __TBB_machine_cmpswp1(P,V,C)
#define __TBB_CompareAndSwap2(P,V,C) __TBB_machine_cmpswp2(P,V,C)
#define __TBB_CompareAndSwap4(P,V,C) __TBB_machine_cmpswp4(P,V,C)
#define __TBB_CompareAndSwap8(P,V,C) __TBB_machine_cmpswp8(P,V,C)
#define __TBB_CompareAndSwapW(P,V,C) __TBB_machine_cmpswp4(P,V,C)

#define __TBB_FetchAndAdd1(P,V) __TBB_machine_fetchadd1(P,V)
#define __TBB_FetchAndAdd2(P,V) __TBB_machine_fetchadd2(P,V)
#define __TBB_FetchAndAdd4(P,V) __TBB_machine_fetchadd4(P,V)
#define __TBB_FetchAndAddW(P,V) __TBB_machine_fetchadd4(P,V)

#define __TBB_FetchAndStore1(P,V) __TBB_machine_fetchstore1(P,V)
#define __TBB_FetchAndStore2(P,V) __TBB_machine_fetchstore2(P,V)
#define __TBB_FetchAndStore4(P,V) __TBB_machine_fetchstore4(P,V)
#define __TBB_FetchAndStoreW(P,V) __TBB_machine_fetchstore4(P,V)

#define __TBB_Store8(P,V) __TBB_machine_store8(P,V)
#define __TBB_Load8(P)    __TBB_machine_load8(P)

#define __TBB_AtomicOR(P,V) __TBB_machine_or(P,V)
#define __TBB_AtomicAND(P,V) __TBB_machine_and(P,V)


// Those we chose not to implement (they will be implemented generically using CMPSWP8)
#undef __TBB_FetchAndAdd8
#undef __TBB_FetchAndStore8

// Definition of other functions
#define __TBB_Pause(V) __TBB_machine_pause(V)
#define __TBB_Log2(V)  __TBB_machine_lg(V)

// Special atomic functions
#define __TBB_FetchAndAddWrelease(P,V) __TBB_FetchAndAddW(P,V)
#define __TBB_FetchAndIncrementWacquire(P) __TBB_FetchAndAddW(P,1)
#define __TBB_FetchAndDecrementWrelease(P) __TBB_FetchAndAddW(P,-1)

// Use generic definitions from tbb_machine.h
#undef __TBB_TryLockByte
#undef __TBB_LockByte
