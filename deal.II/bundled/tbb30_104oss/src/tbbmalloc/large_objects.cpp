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

#include "tbbmalloc_internal.h"

/********* Allocation of large objects ************/


namespace rml {
namespace internal {

static struct LargeBlockCacheStat {
    uintptr_t age;
    size_t cacheSize;
} loCacheStat;

 /*
  * The number of bins to cache large objects.
  */
const uint32_t numLargeBlockBins = 1024; // for 1024 max cached size is near 8MB
 

class CachedBlocksList {
    LargeMemoryBlock *first,
                     *last;
    /* age of an oldest block in the list; equal to last->age, if last defined,
       used for quick cheching it without acquiring the lock. */
    uintptr_t     oldest;
    /* currAge when something was excluded out of list because of the age,
       not because of cache hit */
    uintptr_t     lastCleanedAge;
    /* Current threshold value for the blocks of a particular size. 
       Set on cache miss. */
    intptr_t      ageThreshold;

    MallocMutex   lock;
    /* CachedBlocksList should be placed in zero-initialized memory,
       ctor not needed. */
    CachedBlocksList();
public:
    inline void push(LargeMemoryBlock* ptr);
    inline LargeMemoryBlock* pop();
    void releaseLastIfOld(uintptr_t currAge, size_t size);
};

/*
 * Array of bins with lists of recently freed large objects cached for re-use.
 */
static char globalCachedBlockBinsSpace[sizeof(CachedBlocksList)*numLargeBlockBins];
static CachedBlocksList* globalCachedBlockBins = (CachedBlocksList*)globalCachedBlockBinsSpace;

/*
 * Large Objects are the only objects in the system that begin 
 * on a 16K byte boundary since the blocks used for smaller objects 
 * have the Block structure at each 16K boundary.
 */
static uintptr_t cleanupCacheIfNeed();

void CachedBlocksList::push(LargeMemoryBlock *ptr)
{   
    ptr->prev = NULL;
    ptr->age  = cleanupCacheIfNeed ();

    MallocMutex::scoped_lock scoped_cs(lock);
    ptr->next = first;
    first = ptr;
    if (ptr->next) ptr->next->prev = ptr;
    if (!last) {
        MALLOC_ASSERT(0 == oldest, ASSERT_TEXT);
        oldest = ptr->age;
        last = ptr;
    }
}

LargeMemoryBlock *CachedBlocksList::pop()
{   
    uintptr_t currAge = cleanupCacheIfNeed();
    LargeMemoryBlock *result=NULL;
    {
        MallocMutex::scoped_lock scoped_cs(lock);
        if (first) {
            result = first;
            first = result->next;
            if (first)  
                first->prev = NULL;
            else {
                last = NULL;
                oldest = 0;
            }
        } else {
            /* If cache miss occured, set ageThreshold to twice the difference 
               between current time and last time cache was cleaned. */
            ageThreshold = 2*(currAge - lastCleanedAge);
        }
    }
    return result;
}

void CachedBlocksList::releaseLastIfOld(uintptr_t currAge, size_t size)
{
    LargeMemoryBlock *toRelease = NULL;
 
    /* oldest may be more recent then age, that's why cast to signed type
       was used. age overflow is also processed correctly. */
    if (last && (intptr_t)(currAge - oldest) > ageThreshold) {
        MallocMutex::scoped_lock scoped_cs(lock);
        // double check
        if (last && (intptr_t)(currAge - last->age) > ageThreshold) {
            do {
                last = last->prev;
            } while (last && (intptr_t)(currAge - last->age) > ageThreshold);
            if (last) {
                toRelease = last->next;
                oldest = last->age;
                last->next = NULL;
            } else {
                toRelease = first;
                first = NULL;
                oldest = 0;
            }
            MALLOC_ASSERT( toRelease, ASSERT_TEXT );
            lastCleanedAge = toRelease->age;
        } 
        else 
            return;
    }
    while ( toRelease ) {
        LargeMemoryBlock *helper = toRelease->next;
        removeBackRef(toRelease->backRefIdx);
        freeRawMemory(toRelease, size, toRelease->fromMapMemory);
        toRelease = helper;
    }
}

static uintptr_t cleanupCacheIfNeed ()
{
    /* loCacheStat.age overflow is OK, as we only want difference between 
     * its current value and some recent.
     *
     * Both malloc and free should increment loCacheStat.age, as in 
     * a different case multiple cached blocks would have same age,
     * and accuracy of predictors suffers.
     */
    uintptr_t currAge = (uintptr_t)AtomicIncrement((intptr_t&)loCacheStat.age);

    if ( 0 == currAge % cacheCleanupFreq ) {
        size_t objSize;
        int i;

        for (i = numLargeBlockBins-1, 
             objSize = (numLargeBlockBins-1)*largeBlockCacheStep+blockSize; 
             i >= 0; 
             i--, objSize-=largeBlockCacheStep) {
            /* cached block size on iteration is
             * i*largeBlockCacheStep+blockSize, it seems iterative
             * computation of it improves performance.
             */
            // release from cache blocks that are older than ageThreshold
            globalCachedBlockBins[i].releaseLastIfOld(currAge, objSize);
        }
    }
    return currAge;
}

static LargeMemoryBlock* getCachedLargeBlock (size_t size)
{
    MALLOC_ASSERT( size%largeBlockCacheStep==0, ASSERT_TEXT );
    LargeMemoryBlock *lmb = NULL;
    // blockSize is the minimal alignment and thus the minimal size of a large object.
    size_t idx = (size-minLargeObjectSize)/largeBlockCacheStep;
    if (idx<numLargeBlockBins) {
        lmb = globalCachedBlockBins[idx].pop();
        if (lmb) {
            MALLOC_ITT_SYNC_ACQUIRED(globalCachedBlockBins+idx);
            STAT_increment(getThreadId(), ThreadCommonCounters, allocCachedLargeBlk);
        }
    }
    return lmb;
}

void* mallocLargeObject (size_t size, size_t alignment, bool startupAlloc)
{
    LargeMemoryBlock* lmb;
    size_t headersSize = sizeof(LargeMemoryBlock)+sizeof(LargeObjectHdr);
    size_t allocationSize = alignUp(size+headersSize+alignment, largeBlockCacheStep);

    if (startupAlloc || !(lmb = getCachedLargeBlock(allocationSize))) {
        BackRefIdx backRefIdx;

        if ((backRefIdx = BackRefIdx::newBackRef(/*largeObj=*/true)).isInvalid()) 
            return NULL;
        lmb = (LargeMemoryBlock*)getRawMemory(allocationSize, /*useMapMem=*/startupAlloc);
        if (!lmb) return NULL;
        lmb->fromMapMemory = startupAlloc;
        lmb->backRefIdx = backRefIdx;
        lmb->unalignedSize = allocationSize;
        STAT_increment(getThreadId(), ThreadCommonCounters, allocNewLargeObj);
    }

    void *alignedArea = (void*)alignUp((uintptr_t)lmb+headersSize, alignment);
    LargeObjectHdr *header = (LargeObjectHdr*)alignedArea-1;
    header->memoryBlock = lmb;
    header->backRefIdx = lmb->backRefIdx;
    setBackRef(header->backRefIdx, header);
 
    lmb->objectSize = size;

    MALLOC_ASSERT( isLargeObject(alignedArea), ASSERT_TEXT );
    return alignedArea;
}

static bool freeLargeObjectToCache (LargeMemoryBlock* largeBlock)
{
    size_t size = largeBlock->unalignedSize;
    size_t idx = (size-minLargeObjectSize)/largeBlockCacheStep;
    if (idx<numLargeBlockBins) {
        MALLOC_ASSERT( size%largeBlockCacheStep==0, ASSERT_TEXT );
        MALLOC_ITT_SYNC_RELEASING(globalCachedBlockBins+idx);
        globalCachedBlockBins[idx].push(largeBlock);

        STAT_increment(getThreadId(), ThreadCommonCounters, cacheLargeBlk);
        return true;
    }
    return false;
}

void freeLargeObject (void *object)
{
    LargeObjectHdr *header = (LargeObjectHdr*)object - 1;

    // overwrite backRefIdx to simplify double free detection
    header->backRefIdx = BackRefIdx();
    if (!freeLargeObjectToCache(header->memoryBlock)) {
        removeBackRef(header->memoryBlock->backRefIdx);
        freeRawMemory(header->memoryBlock, header->memoryBlock->unalignedSize, 
                      /*useMapMem=*/ header->memoryBlock->fromMapMemory);
        STAT_increment(getThreadId(), ThreadCommonCounters, freeLargeObj);
    }
}

/*********** End allocation of large objects **********/



} // namespace internal
} // namespace rml

