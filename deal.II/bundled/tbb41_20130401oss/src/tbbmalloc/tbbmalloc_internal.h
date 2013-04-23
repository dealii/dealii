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

#ifndef __TBB_tbbmalloc_internal_H
#define __TBB_tbbmalloc_internal_H 1


#include "TypeDefinitions.h" /* Also includes customization layer Customize.h */

#if USE_PTHREAD
    // Some pthreads documentation says that <pthreads.h> must be first header.
    #include <pthread.h>
    typedef pthread_key_t tls_key_t;
#elif USE_WINTHREAD
    #include "tbb/machine/windows_api.h"
    typedef DWORD tls_key_t;
#else
    #error Must define USE_PTHREAD or USE_WINTHREAD
#endif

#include "tbb/tbb_config.h"
#if __TBB_LIBSTDCPP_EXCEPTION_HEADERS_BROKEN
  #define _EXCEPTION_PTR_H /* prevents exception_ptr.h inclusion */
  #define _GLIBCXX_NESTED_EXCEPTION_H /* prevents nested_exception.h inclusion */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <limits.h> // for CHAR_BIT
#include <string.h> // for memset
#if MALLOC_CHECK_RECURSION
#include <new>        /* for placement new */
#endif
#include "tbb/scalable_allocator.h"
#include "tbbmalloc_internal_api.h"

#if __sun || __SUNPRO_CC
#define __asm__ asm
#endif

/********* Various compile-time options        **************/

#if !__TBB_DEFINE_MIC && __TBB_MIC_NATIVE
 #error Intel(R) Many Integrated Core Compiler does not define __MIC__ anymore.
#endif

#define MALLOC_TRACE 0

#if MALLOC_TRACE
#define TRACEF(x) printf x
#else
#define TRACEF(x) ((void)0)
#endif /* MALLOC_TRACE */

#define ASSERT_TEXT NULL

#define COLLECT_STATISTICS ( MALLOC_DEBUG && MALLOCENV_COLLECT_STATISTICS )
#include "Statistics.h"

/********* End compile-time options        **************/

namespace rml {

namespace internal {

#if __TBB_MALLOC_LOCACHE_STAT
extern intptr_t mallocCalls, cacheHits;
extern intptr_t memAllocKB, memHitKB;
#endif

//! Utility template function to prevent "unused" warnings by various compilers.
template<typename T>
void suppress_unused_warning( const T& ) {}

/********** Various numeric parameters controlling allocations ********/

/*
 * smabSize - the size of a block for allocation of small objects,
 * it must be larger than maxSegregatedObjectSize.
 */
const uintptr_t slabSize = 16*1024;

/*
 * Large blocks cache cleanup frequency.
 * It should be power of 2 for the fast checking.
 */
const unsigned cacheCleanupFreq = 256;

/*
 * Best estimate of cache line size, for the purpose of avoiding false sharing.
 * Too high causes memory overhead, too low causes false-sharing overhead.
 * Because, e.g., 32-bit code might run on a 64-bit system with a larger cache line size,
 * it would probably be better to probe at runtime where possible and/or allow for an environment variable override,
 * but currently this is still used for compile-time layout of class Block, so the change is not entirely trivial.
 */
#if __powerpc64__ || __ppc64__ || __bgp__
const uint32_t estimatedCacheLineSize = 128;
#else
const uint32_t estimatedCacheLineSize =  64;
#endif

/*
 * Alignment of large (>= minLargeObjectSize) objects.
 */
const size_t largeObjectAlignment = estimatedCacheLineSize;

/********** End of numeric parameters controlling allocations *********/

class BlockI;
struct LargeMemoryBlock;
struct ExtMemoryPool;
struct MemRegion;
class FreeBlock;
class TLSData;
class Backend;
class MemoryPool;
extern const uint32_t minLargeObjectSize;

class TLSKey {
    tls_key_t TLS_pointer_key;
public:
    TLSKey();
   ~TLSKey();
    TLSData* getThreadMallocTLS() const;
    void setThreadMallocTLS( TLSData * newvalue );
    TLSData* createTLS(MemoryPool *memPool, Backend *backend);
};

// TODO: make BitMaskBasic more general
// (currenty, it fits BitMaskMin well, but not as suitable for BitMaskMax)
template<unsigned NUM>
class BitMaskBasic {
    static const int SZ = (NUM-1)/(CHAR_BIT*sizeof(uintptr_t))+1;
    static const unsigned WORD_LEN = CHAR_BIT*sizeof(uintptr_t);
    uintptr_t mask[SZ];
protected:
    void set(size_t idx, bool val) {
        MALLOC_ASSERT(idx<NUM, ASSERT_TEXT);

        size_t i = idx / WORD_LEN;
        int pos = WORD_LEN - idx % WORD_LEN - 1;
        if (val)
            AtomicOr(&mask[i], 1ULL << pos);
        else
            AtomicAnd(&mask[i], ~(1ULL << pos));
    }
    int getMinTrue(unsigned startIdx) const {
        size_t idx = startIdx / WORD_LEN;
        uintptr_t curr;
        int pos;

        if (startIdx % WORD_LEN) { // clear bits before startIdx
            pos = WORD_LEN - startIdx % WORD_LEN;
            curr = mask[idx] & ((1ULL<<pos) - 1);
        } else
            curr = mask[idx];

        for (int i=idx; i<SZ; i++, curr=mask[i]) {
            if (-1 != (pos = BitScanRev(curr)))
                return (i+1)*WORD_LEN - pos - 1;
        }
        return -1;
    }
public:
    void reset() { for (int i=0; i<SZ; i++) mask[i] = 0; }
};

template<unsigned NUM>
class BitMaskMin : public BitMaskBasic<NUM> {
public:
    void set(size_t idx, bool val) { BitMaskBasic<NUM>::set(idx, val); }
    int getMinTrue(unsigned startIdx) const {
        return BitMaskBasic<NUM>::getMinTrue(startIdx);
    }
};

template<unsigned NUM>
class BitMaskMax : public BitMaskBasic<NUM> {
public:
    void set(size_t idx, bool val) {
        BitMaskBasic<NUM>::set(NUM - 1 - idx, val);
    }
    int getMaxTrue(unsigned startIdx) const {
        int p = BitMaskBasic<NUM>::getMinTrue(NUM-startIdx-1);
        return -1==p? -1 : (int)NUM - 1 - p;
    }
};

/* cache blocks in range [MinSize; MaxSize) in bins with CacheStep
 TooLargeFactor -- when cache size treated "too large" in comparison to user data size
 OnMissFactor -- If cache miss occurred and cache was cleaned,
                 set ageThreshold to OnMissFactor * the difference
                 between current time and last time cache was cleaned.
 LongWaitFactor -- to detect rarely-used bins and forget about their usage history
*/
template<size_t MIN_SIZE, size_t MAX_SIZE, uint32_t CACHE_STEP, int TOO_LARGE,
         int ON_MISS, int LONG_WAIT>
struct LargeObjectCacheProps {
    static const size_t MinSize = MIN_SIZE, MaxSize = MAX_SIZE;
    static const uint32_t CacheStep = CACHE_STEP;
    static const int TooLargeFactor = TOO_LARGE, OnMissFactor = ON_MISS,
        LongWaitFactor = LONG_WAIT;
};

template<typename Props>
class LargeObjectCacheImpl {
    // The number of bins to cache large objects.
    static const uint32_t numBins = (Props::MaxSize-Props::MinSize)/Props::CacheStep;

    typedef BitMaskMax<numBins> BinBitMask;

    // Current sizes of used and cached objects. It's calculated while we are
    // traversing bins, and used for isLOCTooLarge() check at the same time.
    class BinsSummary {
        size_t usedSz;
        size_t cachedSz;
    public:
        BinsSummary() : usedSz(0), cachedSz(0) {}
        // "too large" criteria
        bool isLOCTooLarge() const { return cachedSz > Props::TooLargeFactor*usedSz; }
        void update(size_t usedSize, size_t cachedSize) {
            usedSz += usedSize;
            cachedSz += cachedSize;
        }
        void reset() { usedSz = cachedSz = 0; }
    };

    // 2-linked list of same-size cached blocks ordered by age (oldest on top)
    // TODO: are we really want the list to be 2-linked? This allows us
    // reduce memory consumption and do less operations under lock.
    // TODO: try to switch to 32-bit logical time to save space in CacheBin
    // and move bins to different cache lines.
    class CacheBin {
        LargeMemoryBlock *first,
                         *last;
  /* age of an oldest block in the list; equal to last->age, if last defined,
     used for quick cheching it without acquiring the lock. */
        uintptr_t         oldest;
  /* currAge when something was excluded out of list because of the age,
     not because of cache hit */
        uintptr_t         lastCleanedAge;
  /* Current threshold value for the blocks of a particular size.
     Set on cache miss. */
        intptr_t          ageThreshold;

  /* total size of all objects corresponding to the bin and allocated by user */
        size_t            usedSize,
  /* total size of all objects cached in the bin */
                          cachedSize;
  /* time of last hit for the bin */
        intptr_t          lastHit;
  /* time of last get called for the bin */
        uintptr_t         lastGet;

        MallocMutex       lock;
  /* should be placed in zero-initialized memory, ctor not needed. */
        CacheBin();
        void forgetOutdatedState(uintptr_t currT);
    public:
        void init() { memset(this, 0, sizeof(CacheBin)); }
        LargeMemoryBlock *putList(ExtMemoryPool *extMemPool, LargeMemoryBlock *head, BinBitMask *bitMask, int idx);
        inline LargeMemoryBlock *get(size_t size, uintptr_t currTime, bool *setNonEmpty);
        void decreaseThreshold() {
            if (ageThreshold)
                ageThreshold = (ageThreshold + lastHit)/2;
        }
        void updateBinsSummary(BinsSummary *binsSummary) const {
            binsSummary->update(usedSize, cachedSize);
        }
        bool cleanToThreshold(Backend *backend, BinBitMask *bitMask, uintptr_t currTime, int idx);
        bool cleanAll(Backend *backend, BinBitMask *bitMask, int idx);
        void decrUsedSize(size_t size, BinBitMask *bitMask, int idx) {
            MallocMutex::scoped_lock scoped_cs(lock);
            usedSize -= size;
            if (!usedSize && !first)
                bitMask->set(idx, false);
        }
        size_t getSize() const { return cachedSize; }
        size_t getUsedSize() const { return usedSize; }
        size_t reportStat(int num, FILE *f);
    };

    intptr_t     tooLargeLOC; // how many times LOC was "too large"
    // for fast finding of used bins and bins with non-zero usedSize;
    // indexed from the end, as we need largest 1st
    BinBitMask   bitMask;
    // bins with lists of recently freed large blocks cached for re-use
    CacheBin bin[numBins];

public:
    static int sizeToIdx(size_t size) {
        MALLOC_ASSERT(Props::MinSize <= size && size < Props::MaxSize, ASSERT_TEXT);
        return (size-Props::MinSize)/Props::CacheStep;
    }
    static int getNumBins() { return numBins; }

    void putList(ExtMemoryPool *extMemPool, LargeMemoryBlock *largeBlock);
    LargeMemoryBlock *get(uintptr_t currTime, size_t size);

    void rollbackCacheState(size_t size);
    uintptr_t cleanupCacheIfNeeded(ExtMemoryPool *extMemPool, uintptr_t currTime);
    bool regularCleanup(Backend *backend, uintptr_t currAge);
    bool cleanAll(Backend *backend);
    void reset() {
        tooLargeLOC = 0;
        for (int i = numBins-1; i >= 0; i--)
            bin[i].init();
        bitMask.reset();
    }
#if __TBB_MALLOC_LOCACHE_STAT
    void reportStat(FILE *f);
#endif
#if __TBB_MALLOC_WHITEBOX_TEST
    size_t getLOCSize() const;
    size_t getUsedSize() const;
#endif
};

class LargeObjectCache {
    static const size_t minLargeSize =  8*1024,
                        maxLargeSize =  8*1024*1024,
                        maxHugeSize = 128*1024*1024;
public:
    // Difference between object sizes in large block bins
    static const uint32_t largeBlockCacheStep =  8*1024,
                          hugeBlockCacheStep = 512*1024;
private:
    typedef LargeObjectCacheImpl< LargeObjectCacheProps<minLargeSize, maxLargeSize, largeBlockCacheStep, 2, 2, 16> > LargeCacheType;
    typedef LargeObjectCacheImpl< LargeObjectCacheProps<maxLargeSize, maxHugeSize, hugeBlockCacheStep, 1, 1, 4> > HugeCacheType;

    LargeCacheType largeCache;
    HugeCacheType hugeCache;

    /* logical time, incremented on each put/get operation
       To prevent starvation between pools, keep separatly for each pool.
       Overflow is OK, as we only want difference between
       its current value and some recent.

       Both malloc and free should increment logical time, as in
       a different case multiple cached blocks would have same age,
       and accuracy of predictors suffers.
    */
    uintptr_t cacheCurrTime;

    static int sizeToIdx(size_t size);
    bool doRegularCleanup(Backend *backend, uintptr_t currTime);
public:
    void put(ExtMemoryPool *extMemPool, LargeMemoryBlock *largeBlock);
    void putList(ExtMemoryPool *extMemPool, LargeMemoryBlock *head);
    LargeMemoryBlock *get(Backend *backend, size_t size);

    void rollbackCacheState(size_t size);
    void cleanupCacheIfNeeded(Backend *backend, uintptr_t currTime);
    void cleanupCacheIfNeededOnRange(Backend *backend, uintptr_t range, uintptr_t currTime);
    bool regularCleanup(Backend *backend) {
        return doRegularCleanup(backend, FencedLoad((intptr_t&)cacheCurrTime));
    }
    bool cleanAll(Backend *backend);
    void reset() {
        largeCache.reset();
        hugeCache.reset();
    }
#if __TBB_MALLOC_LOCACHE_STAT
    void reportStat(FILE *f);
#endif
#if __TBB_MALLOC_WHITEBOX_TEST
    size_t getLOCSize() const;
    size_t getUsedSize() const;
#endif
    static size_t alignToBin(size_t size) {
        return size<maxLargeSize? alignUp(size, largeBlockCacheStep)
            : alignUp(size, hugeBlockCacheStep);
    }

    uintptr_t getCurrTime();
    uintptr_t getCurrTimeRange(uintptr_t range);
};

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

// Block header is used during block coalescing
// and must be preserved in used blocks.
class BlockI {
    intptr_t     blockState[2];
};

struct LargeMemoryBlock : public BlockI {
    LargeMemoryBlock *next,          // ptrs in list of cached blocks
                     *prev,
    // 2-linked list of pool's large objects
    // Used to destroy backrefs on pool destroy/reset (backrefs are global)
    // and for releasing all non-binned blocks.
                     *gPrev,
                     *gNext;
    uintptr_t         age;           // age of block while in cache
    size_t            objectSize;    // the size requested by a client
    size_t            unalignedSize; // the size requested from getMemory
    BackRefIdx        backRefIdx;    // cached here, used copy is in LargeObjectHdr
};

// global state of blocks currently in processing
class BackendSync {
    // Class instances should reside in zero-initialized memory!
    // The number of blocks currently removed from a bin and not returned back
    intptr_t  blocksInProcessing;  // to another
    intptr_t  binsModifications;   // incremented on every bin modification
public:
    void consume() { AtomicIncrement(blocksInProcessing); }
    void pureSignal() { AtomicIncrement(binsModifications); }
    void signal() {
#if __TBB_MALLOC_BACKEND_STAT
        MALLOC_ITT_SYNC_RELEASING(&blocksInProcessing);
#endif
        AtomicIncrement(binsModifications);
        intptr_t prev = AtomicAdd(blocksInProcessing, -1);
        MALLOC_ASSERT(prev > 0, ASSERT_TEXT);
        suppress_unused_warning(prev);
    }
    intptr_t getNumOfMods() const { return FencedLoad(binsModifications); }
    // return true if need re-do the search
    bool waitTillSignalled(intptr_t startModifiedCnt) {
        intptr_t myBlocksNum = FencedLoad(blocksInProcessing);
        if (!myBlocksNum) {
            // no threads, but were bins modified since scanned?
            return startModifiedCnt != getNumOfMods();
        }
#if __TBB_MALLOC_BACKEND_STAT
        MALLOC_ITT_SYNC_PREPARE(&blocksInProcessing);
#endif
        for (;;) {
            SpinWaitWhileEq(blocksInProcessing, myBlocksNum);
            if (myBlocksNum > blocksInProcessing)
                break;
            myBlocksNum = FencedLoad(blocksInProcessing);
        }
#if __TBB_MALLOC_BACKEND_STAT
        MALLOC_ITT_SYNC_ACQUIRED(&blocksInProcessing);
#endif
        return true;
    }
};

class CoalRequestQ { // queue of free blocks that coalescing was delayed
    FreeBlock *blocksToFree;
public:
    FreeBlock *getAll(); // return current list of blocks and make queue empty
    void putBlock(FreeBlock *fBlock);
};

class MemExtendingSema {
    intptr_t     active;
public:
    bool wait() {
        bool rescanBins = false;
        // up to 3 threads can add more memory from OS simultaneously,
        // rest of threads have to wait
        for (;;) {
            intptr_t prevCnt = FencedLoad(active);
            if (prevCnt < 3) {
                intptr_t n = AtomicCompareExchange(active, prevCnt+1, prevCnt);
                if (n == prevCnt)
                    break;
            } else {
                SpinWaitWhileEq(active, prevCnt);
                rescanBins = true;
                break;
            }
        }
        return rescanBins;
    }
    void signal() { AtomicAdd(active, -1); }
};

class Backend {
private:
/* Blocks in range [minBinnedSize; getMaxBinnedSize()] are kept in bins,
   one region can contains several blocks. Larger blocks are allocated directly
   and one region always contains one block.
*/
    enum {
        minBinnedSize = 8*1024UL,
        /*   If huge pages are available, maxBinned_HugePage used.
             If not, maxBinned_SmallPage is the thresold.
             TODO: use pool's granularity for upper bound setting.*/
        maxBinned_SmallPage = 1024*1024UL,
        // TODO: support other page sizes
        maxBinned_HugePage = 4*1024*1024UL
    };
public:
    static const int freeBinsNum =
        (maxBinned_HugePage-minBinnedSize)/LargeObjectCache::largeBlockCacheStep + 1;

    // if previous access missed per-thread slabs pool,
    // allocate numOfSlabAllocOnMiss blocks in advance
    static const int numOfSlabAllocOnMiss = 2;

    enum {
        NO_BIN = -1,
        HUGE_BIN = freeBinsNum-1
    };

    // Bin keeps 2-linked list of free blocks. It must be 2-linked
    // because during coalescing a block it's removed from a middle of the list.
    struct Bin {
        FreeBlock   *head,
                    *tail;
        MallocMutex  tLock;

        void removeBlock(FreeBlock *fBlock);
        void reset() { head = tail = 0; }
#if __TBB_MALLOC_BACKEND_STAT
        size_t countFreeBlocks();
#endif
        bool empty() const { return !head; }
    };

    // array of bins accomplished bitmask for fast finding of non-empty bins
    class IndexedBins {
        BitMaskMin<Backend::freeBinsNum> bitMask;
        Bin                           freeBins[Backend::freeBinsNum];
    public:
        FreeBlock *getBlock(int binIdx, BackendSync *sync, size_t size,
                            bool resSlabAligned, bool alignedBin, bool wait,
                            int *resLocked);
        void lockRemoveBlock(int binIdx, FreeBlock *fBlock);
        void addBlock(int binIdx, FreeBlock *fBlock, size_t blockSz, bool addToTail);
        bool tryAddBlock(int binIdx, FreeBlock *fBlock, bool addToTail);
        int getMinNonemptyBin(unsigned startBin) const {
            int p = bitMask.getMinTrue(startBin);
            return p == -1 ? Backend::freeBinsNum : p;
        }
        void verify();
#if __TBB_MALLOC_BACKEND_STAT
        void reportStat(FILE *f);
#endif
        void reset();
    };
    // number of OS/pool callback calls for more memory
    class AskMemFromOSCounter {
        intptr_t cnt;
    public:
        void OSasked() { AtomicIncrement(cnt); }
        intptr_t get() const { return FencedLoad(cnt); }
    };

private:
    ExtMemoryPool *extMemPool;
    // used for release every region on pool destroying
    MemRegion     *regionList;
    MallocMutex    regionListLock;

    CoalRequestQ   coalescQ; // queue of coalescing requests
    BackendSync    bkndSync;
    // semaphore protecting adding more more memory from OS
    MemExtendingSema memExtendingSema;

    // Using of maximal observed requested size allows descrease
    // memory consumption for small requests and descrease fragmentation
    // for workloads when small and large allocation requests are mixed.
    // TODO: decrease, not only increase it
    size_t         maxRequestedSize;
    void correctMaxRequestSize(size_t requestSize);

    size_t addNewRegion(size_t rawSize, bool exact);
    FreeBlock *findBlockInRegion(MemRegion *region);
    void startUseBlock(MemRegion *region, FreeBlock *fBlock);
    void releaseRegion(MemRegion *region);

    bool askMemFromOS(size_t totalReqSize, intptr_t startModifiedCnt,
                      int *lockedBinsThreshold,
                      int numOfLockedBins, bool *largeBinsUpdated);
    FreeBlock *genericGetBlock(int num, size_t size, bool resSlabAligned);
    void genericPutBlock(FreeBlock *fBlock, size_t blockSz);
    FreeBlock *getFromAlignedSpace(int binIdx, int num, size_t size, bool resSlabAligned, bool wait, int *locked);
    FreeBlock *getFromBin(int binIdx, int num, size_t size, bool resSlabAligned, int *locked);

    FreeBlock *doCoalesc(FreeBlock *fBlock, MemRegion **memRegion);
    void coalescAndPutList(FreeBlock *head, bool forceCoalescQDrop);
    bool scanCoalescQ(bool forceCoalescQDrop);
    void coalescAndPut(FreeBlock *fBlock, size_t blockSz);

    void removeBlockFromBin(FreeBlock *fBlock);

    void *getRawMem(size_t &size) const;
    void freeRawMem(void *object, size_t size) const;

    void putLargeBlock(LargeMemoryBlock *lmb);
public:
    void verify();
#if __TBB_MALLOC_BACKEND_STAT
    void reportStat(FILE *f);
#endif
    bool bootstrap(ExtMemoryPool *extMemoryPool) {
        extMemPool = extMemoryPool;
        return addNewRegion(2*1024*1024, /*exact=*/false);
    }
    void reset();
    bool destroy();

    BlockI *getSlabBlock(int num) {
        BlockI *b = (BlockI*)
            genericGetBlock(num, slabSize, /*resSlabAligned=*/true);
        MALLOC_ASSERT(isAligned(b, slabSize), ASSERT_TEXT);
        return b;
    }
    void putSlabBlock(BlockI *block) {
        genericPutBlock((FreeBlock *)block, slabSize);
    }
    void *getBackRefSpace(size_t size, bool *rawMemUsed);
    void putBackRefSpace(void *b, size_t size, bool rawMemUsed);

    bool inUserPool() const;

    LargeMemoryBlock *getLargeBlock(size_t size);
    void returnLargeObject(LargeMemoryBlock *lmb);

    AskMemFromOSCounter askMemFromOSCounter;
private:
    static int sizeToBin(size_t size) {
        if (size >= maxBinned_HugePage)
            return HUGE_BIN;
        else if (size < minBinnedSize)
            return NO_BIN;

        int bin = (size - minBinnedSize)/LargeObjectCache::largeBlockCacheStep;

        MALLOC_ASSERT(bin < HUGE_BIN, "Invalid size.");
        return bin;
    }
#if __TBB_MALLOC_BACKEND_STAT
    static size_t binToSize(int bin) {
        MALLOC_ASSERT(bin < HUGE_BIN, "Invalid bin.");

        return bin*largeBlockCacheStep + minBinnedSize;
    }
#endif
    static bool toAlignedBin(FreeBlock *block, size_t size) {
        return isAligned((char*)block+size, slabSize)
            && size >= slabSize;
    }
    inline size_t getMaxBinnedSize();

    IndexedBins freeLargeBins,
                freeAlignedBins;
};

class AllLargeBlocksList {
    MallocMutex       largeObjLock;
    LargeMemoryBlock *loHead;
public:
    LargeMemoryBlock *getHead() { return loHead; }
    void add(LargeMemoryBlock *lmb);
    void remove(LargeMemoryBlock *lmb);
    void removeAll(Backend *backend);
};

// An TBB allocator mode that can be controlled by user
// via API/environment variable. Must be placed in zero-initialized memory.
// External synchronization assumed.
// TODO: TBB_VERSION support
class AllocControlledMode {
    intptr_t val;
    bool     setDone;
public:
    intptr_t get() const {
        MALLOC_ASSERT(setDone, ASSERT_TEXT);
        return val;
    }
    void set(intptr_t newVal) { // note set() can be called before init()
        val = newVal;
        setDone = true;
    }
    // envName - environment variable to get controlled mode
    void initReadEnv(const char *envName, intptr_t defaultVal);
};

// init() and printStatus() is called only under global initialization lock.
// Race is possible between registerAllocation() and registerReleasing(),
// harm is that up to single huge page releasing is missed (because failure
// to get huge page is registred only 1st time), that is negligible.
// setMode is also can be called concurrently.
// Object must reside in zero-initialized memory
class HugePagesStatus {
private:
    AllocControlledMode requestedMode; // changed only by user
               // to keep enabled and requestedMode consistent
    MallocMutex setModeLock;
    size_t      pageSize;
    intptr_t    needActualStatusPrint;

    static void doPrintStatus(bool state, const char *stateName);
public:
    // both variables are changed only inside HugePagesStatus
    intptr_t    enabled;
    // Have we got huge pages at all? It's used when large hugepage-aligned
    // region is releasing, to find can it release some huge pages or not.
    intptr_t    wasObserved;

    size_t getSize() const {
        MALLOC_ASSERT(pageSize, ASSERT_TEXT);
        return pageSize;
    }
    void printStatus();
    void registerAllocation(bool available);
    void registerReleasing(size_t size);

    void init(size_t hugePageSize) {
        MALLOC_ASSERT(!hugePageSize || isPowerOfTwo(hugePageSize),
                      "Only memory pages of a power-of-two size are supported.");
        MALLOC_ASSERT(!pageSize, "Huge page size can't be set twice.");
        pageSize = hugePageSize;

        MallocMutex::scoped_lock lock(setModeLock);
        requestedMode.initReadEnv("TBB_MALLOC_USE_HUGE_PAGES", 0);
        enabled = pageSize && requestedMode.get();
    }
    void setMode(intptr_t newVal) {
        MallocMutex::scoped_lock lock(setModeLock);
        requestedMode.set(newVal);
        enabled = pageSize && newVal;
    }
};

extern HugePagesStatus hugePages;

struct ExtMemoryPool {
    Backend           backend;

    intptr_t          poolId;
    // to find all large objects
    AllLargeBlocksList lmbList;
    // Callbacks to be used instead of MapMemory/UnmapMemory.
    rawAllocType      rawAlloc;
    rawFreeType       rawFree;
    size_t            granularity;
    bool              keepAllMemory,
                      delayRegsReleasing,
                      fixedPool;
    TLSKey            tlsPointerKey;  // per-pool TLS key

    LargeObjectCache  loc;

    bool init(intptr_t poolId, rawAllocType rawAlloc, rawFreeType rawFree,
              size_t granularity, bool keepAllMemory, bool fixedPool);
    void initTLS();

    // i.e., not system default pool for scalable_malloc/scalable_free
    bool userPool() const { return rawAlloc; }

     // true if something has beed released
    bool softCachesCleanup();
    bool releaseTLCaches();
    // TODO: to release all thread's pools, not just current thread
    bool hardCachesCleanup();
    void reset() {
        lmbList.removeAll(&backend);
        loc.reset();
        tlsPointerKey.~TLSKey();
        backend.reset();
    }
    void destroy() {
        // pthread_key_dtors must be disabled before memory unmapping
        // TODO: race-free solution
        tlsPointerKey.~TLSKey();
        if (rawFree || !userPool())
            backend.destroy();
    }
    bool mustBeAddedToGlobalLargeBlockList() const { return userPool(); }
    void delayRegionsReleasing(bool mode) { delayRegsReleasing = mode; }
    inline bool regionsAreReleaseable() const;

    LargeMemoryBlock *mallocLargeObject(size_t allocationSize);
    void freeLargeObject(LargeMemoryBlock *lmb);
    void freeLargeObjectList(LargeMemoryBlock *head);
};

inline bool Backend::inUserPool() const { return extMemPool->userPool(); }

struct LargeObjectHdr {
    LargeMemoryBlock *memoryBlock;
    /* Backreference points to LargeObjectHdr.
       Duplicated in LargeMemoryBlock to reuse in subsequent allocations. */
    BackRefIdx       backRefIdx;
};

struct FreeObject {
    FreeObject  *next;
};

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

bool isLargeObject(void *object);

unsigned int getThreadId();

bool initBackRefMaster(Backend *backend);
void destroyBackRefMaster(Backend *backend);
void removeBackRef(BackRefIdx backRefIdx);
void setBackRef(BackRefIdx backRefIdx, void *newPtr);
void *getBackRef(BackRefIdx backRefIdx);

} // namespace internal
} // namespace rml

#endif // __TBB_tbbmalloc_internal_H
