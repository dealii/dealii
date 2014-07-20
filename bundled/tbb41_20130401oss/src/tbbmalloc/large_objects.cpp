/*
    Copyright 2005-2013 Intel Corporation.  All Rights Reserved.

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

#if __TBB_MALLOC_LOCACHE_STAT
intptr_t mallocCalls, cacheHits;
intptr_t memAllocKB, memHitKB;
#endif

inline bool lessThanWithOverflow(intptr_t a, intptr_t b)
{
    return (a < b && (b - a < UINTPTR_MAX/2)) ||
           (a > b && (a - b > UINTPTR_MAX/2));
}

template<typename Props>
LargeMemoryBlock *LargeObjectCacheImpl<Props>::CacheBin::
    putList(ExtMemoryPool *extMemPool, LargeMemoryBlock *head, BinBitMask *bitMask, int idx)
{
    int i, num, totalNum;
    size_t size = head->unalignedSize;
    LargeMemoryBlock *curr, *tail, *toRelease = NULL;
    uintptr_t currTime;

    // we not kept prev pointers during assigning blocks to bins, set them now
    head->prev = NULL;
    for (num=1, curr=head; curr->next; num++, curr=curr->next)
        curr->next->prev = curr;
    tail = curr;
    totalNum = num;

    {
        MallocMutex::scoped_lock scoped_cs(lock);
        usedSize -= num*size;
        // to keep ordering on list, get time under list lock
        currTime = extMemPool->loc.getCurrTimeRange(num);

        for (curr=tail, i=0; curr; curr=curr->prev, i++) {
            curr->age = currTime+i;
            STAT_increment(getThreadId(), ThreadCommonCounters, cacheLargeBlk);
        }

        if (!lastCleanedAge) {
            // 1st object of such size was released.
            // Not cache it, and remeber when this occurs
            // to take into account during cache miss.
            lastCleanedAge = tail->age;
            toRelease = tail;
            tail = tail->prev;
            if (tail)
                tail->next = NULL;
            else
                head = NULL;
            num--;
        }
        if (num) {
            // add [head;tail] list to cache
            tail->next = first;
            if (first)
                first->prev = tail;
            first = head;
            if (!last) {
                MALLOC_ASSERT(0 == oldest, ASSERT_TEXT);
                oldest = tail->age;
                last = tail;
            }

            cachedSize += num*size;
        }
/* It's accebtable, if a bin is empty, and we have non-empty in bit mask.
   So set true in bitmask without lock.
   It's not acceptable, if a bin is non-empty and we have empty in bitmask.
   So set false in bitmask under lock. */

        // No used object, and nothing in the bin, mark the bin as empty
        if (!usedSize && !first)
            bitMask->set(idx, false);
    }
    extMemPool->loc.cleanupCacheIfNeededOnRange(&extMemPool->backend, totalNum, currTime);
    if (toRelease)
        toRelease->prev = toRelease->next = NULL;
    return toRelease;
}

template<typename Props>
LargeMemoryBlock *LargeObjectCacheImpl<Props>::CacheBin::
    get(size_t size, uintptr_t currTime, bool *setNonEmpty)
{
    LargeMemoryBlock *result=NULL;
    {
        MallocMutex::scoped_lock scoped_cs(lock);
        forgetOutdatedState(currTime);

        if (first) {
            result = first;
            first = result->next;
            if (first)
                first->prev = NULL;
            else {
                last = NULL;
                oldest = 0;
            }
            // use moving average with current hit interval
            intptr_t hitR = currTime - result->age;
            lastHit = lastHit? (lastHit + hitR)/2 : hitR;

            cachedSize -= size;
        } else {
            if (lastCleanedAge)
                ageThreshold = Props::OnMissFactor*(currTime - lastCleanedAge);
        }
        if (!usedSize) // inform that there are used blocks in the bin
            *setNonEmpty = true;
        // subject to later correction, if got cache miss and later allocation failed
        usedSize += size;
        lastGet = currTime;
    }
    return result;
}

// forget the history for the bin if it was unused for long time
template<typename Props>
void LargeObjectCacheImpl<Props>::CacheBin::forgetOutdatedState(uintptr_t currTime)
{
    // If the time since the last get is LongWaitFactor times more than ageThreshold
    // for the bin, treat the bin as rarely-used and forget everything we know
    // about it.
    // If LongWaitFactor is too small, we forget too early and
    // so prevents good caching, while if too high, caching blocks
    // with unrelated usage pattern occurs.
    const uintptr_t sinceLastGet = currTime - lastGet;
    bool doCleanup = false;

    if (!last) { // clean only empty bins
        if (ageThreshold)
            doCleanup = sinceLastGet > Props::LongWaitFactor*ageThreshold;
        else if (lastCleanedAge)
            doCleanup = sinceLastGet > Props::LongWaitFactor*(lastCleanedAge - lastGet);
    }
    if (doCleanup) {
        lastCleanedAge = 0;
        ageThreshold = 0;
    }
}

template<typename Props>
bool LargeObjectCacheImpl<Props>::CacheBin::
    cleanToThreshold(Backend *backend, BinBitMask *bitMask, uintptr_t currTime, int idx)
{
    LargeMemoryBlock *toRelease = NULL;
    bool released = false;
#if MALLOC_DEBUG
    uintptr_t nextAge = 0;
#endif

    /* oldest may be more recent then age, that's why cast to signed type
       was used. age overflow is also processed correctly. */
    if (last && (intptr_t)(currTime - oldest) > ageThreshold) {
        MallocMutex::scoped_lock scoped_cs(lock);
        // double check
        if (last && (intptr_t)(currTime - last->age) > ageThreshold) {
            do {
#if MALLOC_DEBUG
                // check that list ordered
                MALLOC_ASSERT(!nextAge || lessThanWithOverflow(nextAge, last->age),
                              ASSERT_TEXT);
                nextAge = last->age;
#endif
                cachedSize -= last->unalignedSize;
                last = last->prev;
            } while (last && (intptr_t)(currTime - last->age) > ageThreshold);
            if (last) {
                toRelease = last->next;
                oldest = last->age;
                last->next = NULL;
            } else {
                toRelease = first;
                first = NULL;
                oldest = 0;
                if (!usedSize)
                    bitMask->set(idx, false);
            }
            MALLOC_ASSERT( toRelease, ASSERT_TEXT );
            lastCleanedAge = toRelease->age;
        }
        else
            return false;
    }
    released = toRelease;

    while ( toRelease ) {
        LargeMemoryBlock *helper = toRelease->next;
        backend->returnLargeObject(toRelease);
        toRelease = helper;
    }
    return released;
}

template<typename Props>
bool LargeObjectCacheImpl<Props>::
    CacheBin::cleanAll(Backend *backend, BinBitMask *bitMask, int idx)
{
    LargeMemoryBlock *toRelease = NULL;
    bool released = false;

    if (last) {
        MallocMutex::scoped_lock scoped_cs(lock);
        // double check
        if (last) {
            toRelease = first;
            last = NULL;
            first = NULL;
            oldest = 0;
            cachedSize = 0;
            if (!usedSize)
                bitMask->set(idx, false);
        }
        else
            return false;
    }
    released = toRelease;

    while ( toRelease ) {
        LargeMemoryBlock *helper = toRelease->next;
        MALLOC_ASSERT(!helper || lessThanWithOverflow(helper->age, toRelease->age),
                      ASSERT_TEXT);
        backend->returnLargeObject(toRelease);
        toRelease = helper;
    }
    return released;
}

template<typename Props>
size_t LargeObjectCacheImpl<Props>::CacheBin::reportStat(int num, FILE *f)
{
#if __TBB_MALLOC_LOCACHE_STAT
    if (first)
        printf("%d(%lu): total %lu KB thr %ld lastCln %lu lastHit %lu oldest %lu\n",
               num, num*CacheStep+MinSize,
               cachedSize/1024, ageThreshold, lastCleanedAge, lastHit, oldest);
#else
    suppress_unused_warning(num);
    suppress_unused_warning(f);
#endif
    return cachedSize;
}

// release from cache blocks that are older than ageThreshold
template<typename Props>
bool LargeObjectCacheImpl<Props>::regularCleanup(Backend *backend, uintptr_t currTime)
{
    bool released = false, doThreshDecr = false;
    BinsSummary binsSummary;

    for (int i = bitMask.getMaxTrue(numBins-1); i >= 0;
         i = bitMask.getMaxTrue(i-1)) {
        bin[i].updateBinsSummary(&binsSummary);
        if (!doThreshDecr && tooLargeLOC>2 && binsSummary.isLOCTooLarge()) {
            // if LOC is too large for quite long time, decrease the threshold
            // based on bin hit statistics.
            // For this, redo cleanup from the beginnig.
            // Note: on this iteration total usedSz can be not too large
            // in comparison to total cachedSz, as we calculated it only
            // partially. We are ok this it.
            i = bitMask.getMaxTrue(numBins-1);
            doThreshDecr = true;
            binsSummary.reset();
            continue;
        }
        if (doThreshDecr)
            bin[i].decreaseThreshold();
        if (bin[i].cleanToThreshold(backend, &bitMask, currTime, i))
            released = true;
    }

    // We want to find if LOC was too large for some time continuously,
    // so OK with races between incrementing and zeroing, but incrementing
    // must be atomic.
    if (binsSummary.isLOCTooLarge())
        AtomicIncrement(tooLargeLOC);
    else
        tooLargeLOC = 0;
    return released;
}

template<typename Props>
bool LargeObjectCacheImpl<Props>::cleanAll(Backend *backend)
{
    bool released = false;
    for (int i = numBins-1; i >= 0; i--)
        released |= bin[i].cleanAll(backend, &bitMask, i);
    return released;
}

#if __TBB_MALLOC_WHITEBOX_TEST
template<typename Props>
size_t LargeObjectCacheImpl<Props>::getLOCSize() const
{
    size_t size = 0;
    for (int i = numBins-1; i >= 0; i--)
        size += bin[i].getSize();
    return size;
}

size_t LargeObjectCache::getLOCSize() const
{
    return largeCache.getLOCSize() + hugeCache.getLOCSize();
}

template<typename Props>
size_t LargeObjectCacheImpl<Props>::getUsedSize() const
{
    size_t size = 0;
    for (int i = numBins-1; i >= 0; i--)
        size += bin[i].getUsedSize();
    return size;
}

size_t LargeObjectCache::getUsedSize() const
{
    return largeCache.getUsedSize() + hugeCache.getUsedSize();
}
#endif // __TBB_MALLOC_WHITEBOX_TEST

uintptr_t LargeObjectCache::getCurrTime()
{
    return (uintptr_t)AtomicIncrement((intptr_t&)cacheCurrTime);
}

uintptr_t LargeObjectCache::getCurrTimeRange(uintptr_t range)
{
    return (uintptr_t)AtomicAdd((intptr_t&)cacheCurrTime, range)+1;
}

void LargeObjectCache::cleanupCacheIfNeeded(Backend *backend, uintptr_t currTime)
{
    if ( 0 == currTime % cacheCleanupFreq )
        doRegularCleanup(backend, currTime);
}

void LargeObjectCache::
    cleanupCacheIfNeededOnRange(Backend *backend, uintptr_t range, uintptr_t currTime)
{
    if (range >= cacheCleanupFreq
        || currTime+range < currTime-1 // overflow, 0 is power of 2, do cleanup
        // (prev;prev+range] contains n*cacheCleanupFreq
        || alignUp(currTime, cacheCleanupFreq)<=currTime+range)
        doRegularCleanup(backend, currTime);
}

bool LargeObjectCache::doRegularCleanup(Backend *backend, uintptr_t currTime)
{
    return largeCache.regularCleanup(backend, currTime)
        | hugeCache.regularCleanup(backend, currTime);
}

bool LargeObjectCache::cleanAll(Backend *backend)
{
    return largeCache.cleanAll(backend) | hugeCache.cleanAll(backend);
}

template<typename Props>
LargeMemoryBlock *LargeObjectCacheImpl<Props>::get(uintptr_t currTime, size_t size)
{
    MALLOC_ASSERT( size%Props::CacheStep==0, ASSERT_TEXT );
    int idx = sizeToIdx(size);
    bool setNonEmpty = false;

    LargeMemoryBlock *lmb = bin[idx].get(size, currTime, &setNonEmpty);
    // Setting to true is possible out of lock. As bitmask is used only for cleanup,
    // the lack of consistency is not violating correctness here.
    if (setNonEmpty)
        bitMask.set(idx, true);
    if (lmb) {
        MALLOC_ITT_SYNC_ACQUIRED(bin+idx);
        STAT_increment(getThreadId(), ThreadCommonCounters, allocCachedLargeBlk);
    }
    return lmb;
}

template<typename Props>
void LargeObjectCacheImpl<Props>::rollbackCacheState(size_t size)
{
    int idx = sizeToIdx(size);
    MALLOC_ASSERT(idx<numBins, ASSERT_TEXT);
    bin[idx].decrUsedSize(size, &bitMask, idx);
}

#if __TBB_MALLOC_LOCACHE_STAT
template<typename Props>
void LargeObjectCacheImpl<Props>::reportStat(FILE *f)
{
    size_t cachedSize = 0;
    for (int i=0; i<numLargeBlockBins; i++)
        cachedSize += bin[i].reportStat(i, f);
    fprintf(f, "total LOC size %lu MB\nnow %lu\n", cachedSize/1024/1024,
            loCacheStat.age);
}

void LargeObjectCache::reportStat(FILE *f)
{
    largeObjs.reportStat(f);
    hugeObjs.reportStat(f);
}
#endif

template<typename Props>
void LargeObjectCacheImpl<Props>::putList(ExtMemoryPool *extMemPool, LargeMemoryBlock *toCache)
{
    int toBinIdx = sizeToIdx(toCache->unalignedSize);

    MALLOC_ITT_SYNC_RELEASING(bin+toBinIdx);
    if (LargeMemoryBlock *release = bin[toBinIdx].putList(extMemPool, toCache,
                                                          &bitMask, toBinIdx))
        extMemPool->backend.returnLargeObject(release);
}

void LargeObjectCache::rollbackCacheState(size_t size)
{
    if (size < maxLargeSize)
        largeCache.rollbackCacheState(size);
    else if (size < maxHugeSize)
        hugeCache.rollbackCacheState(size);
}

// return artifical bin index, it's used only during sorting and never saved
int LargeObjectCache::sizeToIdx(size_t size)
{
    MALLOC_ASSERT(size < maxHugeSize, ASSERT_TEXT);
    return size < maxLargeSize?
        LargeCacheType::sizeToIdx(size) :
        LargeCacheType::getNumBins()+HugeCacheType::sizeToIdx(size);
}

void LargeObjectCache::putList(ExtMemoryPool *extMemPool, LargeMemoryBlock *list)
{
    LargeMemoryBlock *toProcess, *n;

    for (LargeMemoryBlock *curr = list; curr; curr = toProcess) {
        LargeMemoryBlock *tail = curr;
        toProcess = curr->next;
        if (curr->unalignedSize >= maxHugeSize) {
            extMemPool->backend.returnLargeObject(curr);
            continue;
        }
        int currIdx = sizeToIdx(curr->unalignedSize);

        // Find all blocks fitting to same bin. Not use more efficient sorting
        // algorithm because list is short (commonly,
        // LocalLOC's HIGH_MARK-LOW_MARK, i.e. 24 items).
        for (LargeMemoryBlock *b = toProcess; b; b = n) {
            n = b->next;
            if (sizeToIdx(b->unalignedSize) == currIdx) {
                tail->next = b;
                tail = b;
                if (toProcess == b)
                    toProcess = toProcess->next;
                else {
                    b->prev->next = b->next;
                    if (b->next)
                        b->next->prev = b->prev;
                }
            }
        }
        tail->next = NULL;
        if (curr->unalignedSize < maxLargeSize)
            largeCache.putList(extMemPool, curr);
        else
            hugeCache.putList(extMemPool, curr);
    }
}

void LargeObjectCache::put(ExtMemoryPool *extMemPool, LargeMemoryBlock *largeBlock)
{
    if (largeBlock->unalignedSize < maxHugeSize) {
        largeBlock->next = NULL;
        if (largeBlock->unalignedSize<maxLargeSize)
            largeCache.putList(extMemPool, largeBlock);
        else
            hugeCache.putList(extMemPool, largeBlock);
    } else
        extMemPool->backend.returnLargeObject(largeBlock);
}

LargeMemoryBlock *LargeObjectCache::get(Backend *backend, size_t size)
{
    MALLOC_ASSERT( size%largeBlockCacheStep==0, ASSERT_TEXT );
    MALLOC_ASSERT( size>=minLargeSize, ASSERT_TEXT );

    if ( size < maxHugeSize) {
        uintptr_t currTime = getCurrTime();
        cleanupCacheIfNeeded(backend, currTime);
        return size < maxLargeSize?
            largeCache.get(currTime, size) : hugeCache.get(currTime, size);
    }
    return NULL;
}


LargeMemoryBlock *ExtMemoryPool::mallocLargeObject(size_t allocationSize)
{
#if __TBB_MALLOC_LOCACHE_STAT
    AtomicIncrement(mallocCalls);
    AtomicAdd(memAllocKB, allocationSize/1024);
#endif
    LargeMemoryBlock* lmb = loc.get(&backend, allocationSize);
    if (!lmb) {
        BackRefIdx backRefIdx = BackRefIdx::newBackRef(/*largeObj=*/true);
        if (backRefIdx.isInvalid())
            return NULL;

        // unalignedSize is set in getLargeBlock
        lmb = backend.getLargeBlock(allocationSize);
        if (!lmb) {
            removeBackRef(backRefIdx);
            loc.rollbackCacheState(allocationSize);
            return NULL;
        }
        lmb->backRefIdx = backRefIdx;
        STAT_increment(getThreadId(), ThreadCommonCounters, allocNewLargeObj);
    } else {
#if __TBB_MALLOC_LOCACHE_STAT
        AtomicIncrement(cacheHits);
        AtomicAdd(memHitKB, allocationSize/1024);
#endif
    }
    return lmb;
}

void ExtMemoryPool::freeLargeObject(LargeMemoryBlock *mBlock)
{
    loc.put(this, mBlock);
}

void ExtMemoryPool::freeLargeObjectList(LargeMemoryBlock *head)
{
    loc.putList(this, head);
}

bool ExtMemoryPool::softCachesCleanup()
{
    // TODO: cleanup small objects as well
    return loc.regularCleanup(&backend);
}

bool ExtMemoryPool::hardCachesCleanup()
{
    // thread-local caches must be cleaned before LOC,
    // because object from thread-local cache can be released to LOC
    bool tlCaches = releaseTLCaches(), locCaches = loc.cleanAll(&backend);
    return tlCaches || locCaches;
}


/*********** End allocation of large objects **********/

} // namespace internal
} // namespace rml

