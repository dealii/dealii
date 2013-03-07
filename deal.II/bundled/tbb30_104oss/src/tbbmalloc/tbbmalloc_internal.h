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

#ifndef __TBB_tbbmalloc_internal_H
#define __TBB_tbbmalloc_internal_H 1


#include "TypeDefinitions.h" /* Also includes customization layer Customize.h */

#if USE_PTHREAD
    // Some pthreads documentation says that <pthreads.h> must be first header.
    #include <pthread.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#if MALLOC_CHECK_RECURSION
#include <new>        /* for placement new */
#endif

#if __sun || __SUNPRO_CC
#define __asm__ asm 
#endif

extern "C" {
    void * scalable_malloc(size_t size);
    void   scalable_free(void *object);
    void mallocThreadShutdownNotification(void*);
}


/********* Various compile-time options        **************/

#define MALLOC_TRACE 0

#if MALLOC_TRACE
#define TRACEF(x) printf x
#else
#define TRACEF(x) ((void)0)
#endif /* MALLOC_TRACE */

#define ASSERT_TEXT NULL

#define COLLECT_STATISTICS MALLOC_DEBUG && defined(MALLOCENV_COLLECT_STATISTICS)
#include "Statistics.h"

/********* End compile-time options        **************/

namespace rml {
namespace internal {

/********** Various numeric parameters controlling allocations ********/

/*
 * blockSize - the size of a block, it must be larger than maxSegregatedObjectSize.
 *
 */
const uintptr_t blockSize = 16*1024;

/*
 * Difference between object sizes in large block bins
 */
const uint32_t largeBlockCacheStep = 8*1024;

/*
 * Large blocks cache cleanup frequency.
 * It should be power of 2 for the fast checking.
 */
const unsigned cacheCleanupFreq = 256;

/*
 * Alignment of large (>= minLargeObjectSize) objects.
 */
static int largeObjectAlignment = 64; // 64 is common cache line size

/********** End of numeric parameters controlling allocations *********/

class BackRefIdx { // composite index to backreference array
private:
    uint16_t master;      // index in BackRefMaster
    uint16_t largeObj:1;  // is this object "large"?
    uint16_t offset  :15; // offset from beginning of BackRefBlock
public:
    BackRefIdx() : master((uint16_t)-1) {}
    bool isInvalid() const { return master == (uint16_t)-1; }
    bool isLargeObject() const { return largeObj; }
    uint16_t getMaster() const { return master; }
    uint16_t getOffset() const { return offset; }

    // only newBackRef can modify BackRefIdx
    static BackRefIdx newBackRef(bool largeObj);
};

struct LargeMemoryBlock {
    LargeMemoryBlock *next,          // ptrs in list of cached blocks
                     *prev;
    uintptr_t         age;           // age of block while in cache
    size_t            objectSize;    // the size requested by a client
    size_t            unalignedSize; // the size requested from getMemory
    bool              fromMapMemory;
    BackRefIdx        backRefIdx;    // cached here, used copy is in LargeObjectHdr
};

struct LargeObjectHdr {
    LargeMemoryBlock *memoryBlock;
    /* Backreference points to LargeObjectHdr. 
       Duplicated in LargeMemoryBlock to reuse in subsequent allocations. */
    BackRefIdx       backRefIdx;
};

struct FreeObject {
    FreeObject  *next;
};

// interface class for external access to Block
class BlockI {
public:
    static BlockI *getRawBlock(bool startup);
    void initialize(void *bumpPtr);
};

class FreeBlocks {
    typedef void* (*RawAlloc) (size_t size, bool useMapMem);
    typedef void (*RawFree) (void *object, size_t size, bool useMapMem);

    RawAlloc rawAlloc;
    RawFree rawFree;
    size_t memReqSize;

    bool mallocBigBlock();
public:
    bool bootstrap(RawAlloc myAlloc, RawFree myFree, size_t myReqSize);
    BlockI *get(bool startup);
    void put(BlockI *block, bool startup);
    void putList(BlockI *head, BlockI *tail);
};

extern FreeBlocks freeBlocks;

/******* A helper class to support overriding malloc with scalable_malloc *******/
#if MALLOC_CHECK_RECURSION

class RecursiveMallocCallProtector {
    // pointer to an automatic data of holding thread
    static void       *autoObjPtr;
    static MallocMutex rmc_mutex;
    static pthread_t   owner_thread;
/* Under FreeBSD 8.0 1st call to any pthread function including pthread_self
   leads to pthread initialization, that causes malloc calls. As 1st usage of
   RecursiveMallocCallProtector can be before pthread initialized, pthread calls
   can't be used in 1st instance of RecursiveMallocCallProtector.
   RecursiveMallocCallProtector is used 1st time in checkInitialization(),
   so there is a guarantee that on 2nd usage pthread is initialized. 
   No such situation observed with other supported OSes.
 */
#if __FreeBSD__
    static bool        canUsePthread;
#else
    static const bool  canUsePthread = true;
#endif
/*
  The variable modified in checkInitialization,
  so can be read without memory barriers.
 */
    static bool mallocRecursionDetected;

    MallocMutex::scoped_lock* lock_acquired;
    char scoped_lock_space[sizeof(MallocMutex::scoped_lock)+1];

    static uintptr_t absDiffPtr(void *x, void *y) {
        uintptr_t xi = (uintptr_t)x, yi = (uintptr_t)y;
        return xi > yi ? xi - yi : yi - xi;
    }
public:

    RecursiveMallocCallProtector() : lock_acquired(NULL) {
        lock_acquired = new (scoped_lock_space) MallocMutex::scoped_lock( rmc_mutex );
        if (canUsePthread)
            owner_thread = pthread_self();
        autoObjPtr = &scoped_lock_space;
    }
    ~RecursiveMallocCallProtector() {
        if (lock_acquired) {
            autoObjPtr = NULL;
            lock_acquired->~scoped_lock();
        }
    }
    static bool sameThreadActive() {
        if (!autoObjPtr) // fast path
            return false;
        // Some thread has an active recursive call protector; check if the current one.
        // Exact pthread_self based test
        if (canUsePthread) {
            if (pthread_equal( owner_thread, pthread_self() )) {
                mallocRecursionDetected = true;
                return true;
            } else
                return false;
        }
        // inexact stack size based test
        const uintptr_t threadStackSz = 2*1024*1024;
        int dummy;
        return absDiffPtr(autoObjPtr, &dummy)<threadStackSz;
    }
    static bool noRecursion();
/* The function is called on 1st scalable_malloc call to check if malloc calls
   scalable_malloc (nested call must set mallocRecursionDetected). */
    static void detectNaiveOverload() {
        if (!malloc_proxy) {
#if __FreeBSD__
/* If !canUsePthread, we can't call pthread_self() before, but now pthread 
   is already on, so can do it. False positives here lead to silent switching 
   from malloc to mmap for all large allocations with bad performance impact. */
            if (!canUsePthread) {
                canUsePthread = true;
                owner_thread = pthread_self();
            }
#endif
            free(malloc(1));
        }
    }
};

#else

class RecursiveMallocCallProtector {
public:
    RecursiveMallocCallProtector() {}
    ~RecursiveMallocCallProtector() {}
};

#endif  /* MALLOC_CHECK_RECURSION */

bool isMallocInitializedExt();

void* getRawMemory (size_t size, bool useMapMem);
void freeRawMemory (void *object, size_t size, bool useMapMem);

extern const uint32_t minLargeObjectSize;
bool isLargeObject(void *object);
void* mallocLargeObject (size_t size, size_t alignment, bool startupAlloc = false);
void freeLargeObject (void *object);

unsigned int getThreadId();

bool initBackRefMaster();
void removeBackRef(BackRefIdx backRefIdx);
void setBackRef(BackRefIdx backRefIdx, void *newPtr);
void *getBackRef(BackRefIdx backRefIdx);

} // namespace internal
} // namespace rml

#endif // __TBB_tbbmalloc_internal_H
