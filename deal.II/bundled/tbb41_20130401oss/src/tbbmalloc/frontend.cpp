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
#include <errno.h>
#include <new>        /* for placement new */
#include <string.h>   /* for memset */

//! Define the main synchronization method
#define FINE_GRAIN_LOCKS

#if USE_PTHREAD
    #define TlsSetValue_func pthread_setspecific
    #define TlsGetValue_func pthread_getspecific
    #include <sched.h>
    inline void do_yield() {sched_yield();}
    extern "C" { static void mallocThreadShutdownNotification(void*); }

#elif USE_WINTHREAD
#if __TBB_WIN8UI_SUPPORT
#include<thread>
    #define TlsSetValue_func FlsSetValue
    #define TlsGetValue_func FlsGetValue
    #define TlsAlloc() FlsAlloc(NULL)
    #define TlsFree FlsFree
    inline void do_yield() {std::this_thread::yield();}
#else
    #define TlsSetValue_func TlsSetValue
    #define TlsGetValue_func TlsGetValue
    inline void do_yield() {SwitchToThread();}
#endif
#else
    #error Must define USE_PTHREAD or USE_WINTHREAD

#endif


#define FREELIST_NONBLOCKING 1

namespace rml {
class MemoryPool;
namespace internal {

class Block;
class MemoryPool;

#if MALLOC_CHECK_RECURSION

inline bool isMallocInitialized();

bool RecursiveMallocCallProtector::noRecursion() {
    MALLOC_ASSERT(isMallocInitialized(),
                  "Recursion status can be checked only when initialization was done.");
    return !mallocRecursionDetected;
}

#endif // MALLOC_CHECK_RECURSION

/*
 * Block::objectSize value used to mark blocks allocated by startupAlloc
 */
const uint16_t startupAllocObjSizeMark = ~(uint16_t)0;

/*
 * This number of bins in the TLS that leads to blocks that we can allocate in.
 */
const uint32_t numBlockBinLimit = 31;

/*
 * The following constant is used to define the size of struct Block, the block header.
 * The intent is to have the size of a Block multiple of the cache line size, this allows us to
 * get good alignment at the cost of some overhead equal to the amount of padding included in the Block.
 */
const int blockHeaderAlignment = estimatedCacheLineSize;

/********* The data structures and global objects        **************/

/*
 * The malloc routines themselves need to be able to occasionally malloc some space,
 * in order to set up the structures used by the thread local structures. This
 * routine preforms that fuctions.
 */
class BootStrapBlocks {
    MallocMutex bootStrapLock;
    Block      *bootStrapBlock;
    Block      *bootStrapBlockUsed;
    FreeObject *bootStrapObjectList;
public:
    void *allocate(MemoryPool *memPool, size_t size);
    void free(void* ptr);
    void reset();
};

class ThreadId {
    static tls_key_t Tid_key;
    static intptr_t ThreadIdCount;

    unsigned int id;
public:

    static void init() {
#if USE_WINTHREAD
        Tid_key = TlsAlloc();
#else
        int status = pthread_key_create( &Tid_key, NULL );
        if ( status ) {
            fprintf (stderr, "The memory manager cannot create tls key during initialization; exiting \n");
            exit(1);
        }
#endif /* USE_WINTHREAD */
    }
    static void destroy() {
        if( Tid_key ) {
#if USE_WINTHREAD
            TlsFree( Tid_key );
#else
            int status = pthread_key_delete( Tid_key );
            if ( status ) {
                fprintf (stderr, "The memory manager cannot delete tls key; exiting \n");
                exit(1);
            }
#endif /* USE_WINTHREAD */
            Tid_key = 0;
        }
    }
    static ThreadId get() {
        ThreadId result;
        result.id = reinterpret_cast<intptr_t>(TlsGetValue_func(Tid_key));
        if( !result.id ) {
            RecursiveMallocCallProtector scoped;
            // Thread-local value is zero -> first call from this thread,
            // need to initialize with next ID value (IDs start from 1)
            result.id = AtomicIncrement(ThreadIdCount); // returned new value!
            TlsSetValue_func( Tid_key, reinterpret_cast<void*>(result.id) );
        }
        return result;
    }

    bool defined() const { return id; }
    void undef() { id = 0; }
    void invalid() { id = (unsigned int)-1; }
    bool own() const { return id == ThreadId::get().id; }

    friend bool operator==(const ThreadId &id1, const ThreadId &id2);
    friend unsigned int getThreadId();
};

tls_key_t ThreadId::Tid_key;
intptr_t ThreadId::ThreadIdCount;

bool operator==(const ThreadId &id1, const ThreadId &id2) {
    return id1.id == id2.id;
}

unsigned int getThreadId() { return ThreadId::get().id; }

/*********** Code to provide thread ID and a thread-local void pointer **********/

TLSKey::TLSKey()
{
#if USE_WINTHREAD
    TLS_pointer_key = TlsAlloc();
#else
    int status = pthread_key_create( &TLS_pointer_key, mallocThreadShutdownNotification );
    if ( status ) {
        fprintf (stderr, "The memory manager cannot create tls key during initialization; exiting \n");
        exit(1);
    }
#endif /* USE_WINTHREAD */
}

TLSKey::~TLSKey()
{
#if USE_WINTHREAD
    TlsFree(TLS_pointer_key);
#else
    int status1 = pthread_key_delete(TLS_pointer_key);
    if ( status1 ) {
        fprintf (stderr, "The memory manager cannot delete tls key during; exiting \n");
        exit(1);
    }
#endif /* USE_WINTHREAD */
}

inline TLSData* TLSKey::getThreadMallocTLS() const
{
    return (TLSData *)TlsGetValue_func( TLS_pointer_key );
}

inline void TLSKey::setThreadMallocTLS( TLSData * newvalue ) {
    RecursiveMallocCallProtector scoped;
    TlsSetValue_func( TLS_pointer_key, newvalue );
}

/* The 'next' field in the block header has to maintain some invariants:
 *   it needs to be on a 16K boundary and the first field in the block.
 *   Any value stored there needs to have the lower 14 bits set to 0
 *   so that various assert work. This means that if you want to smash this memory
 *   for debugging purposes you will need to obey this invariant.
 * The total size of the header needs to be a power of 2 to simplify
 * the alignment requirements. For now it is a 128 byte structure.
 * To avoid false sharing, the fields changed only locally are separated
 * from the fields changed by foreign threads.
 * Changing the size of the block header would require to change
 * some bin allocation sizes, in particular "fitting" sizes (see above).
 */
class Bin;
class StartupBlock;
class TLSData;

class LifoList {
public:
    inline LifoList();
    inline void push(Block *block);
    inline Block *pop();

private:
    Block *top;
#ifdef FINE_GRAIN_LOCKS
    MallocMutex lock;
#endif /* FINE_GRAIN_LOCKS     */
};

/*
 * When a block that is not completely free is returned for reuse by other threads
 * this is where the block goes.
 *
 * LifoList assumes zero initialization; so below its constructors are omitted,
 * to avoid linking with C++ libraries on Linux.
 */

class OrphanedBlocks {
    LifoList bins[numBlockBinLimit];
public:
    Block *get(Bin* bin, unsigned int size);
    void put(Bin* bin, Block *block);
    void reset();
};

class MemoryPool {
    // if no explicit grainsize, expect to see malloc in user's pAlloc
    // and set reasonable low granularity
    static const size_t defaultGranularity = estimatedCacheLineSize;

    MemoryPool();                  // deny
public:
    static MallocMutex  memPoolListLock;

    // list of all active pools is used to release
    // all TLS data on thread termination or library unload
    MemoryPool    *next,
                  *prev;
    ExtMemoryPool  extMemPool;
    OrphanedBlocks orphanedBlocks;
    BootStrapBlocks bootStrapBlocks;

    bool init(intptr_t poolId, const MemPoolPolicy* memPoolPolicy);
    static void initDefaultPool();
    void reset();
    void destroy();
    void processThreadShutdown(TLSData *tlsData);

    inline TLSData *getTLS(bool create);
    void clearTLS() { extMemPool.tlsPointerKey.setThreadMallocTLS(NULL); }

    Bin *getAllocationBin(TLSData* tls, size_t size);
    Block *getEmptyBlock(size_t size);
    void returnEmptyBlock(Block *block, bool poolTheBlock);

    // get/put large object to/from local large object cache
    void *getFromLLOCache(TLSData *tls, size_t size, size_t alignment);
    void putToLLOCache(TLSData *tls, void *object);

    inline void allocatorCalledHook(TLSData *tls);
};

static char defaultMemPool_space[sizeof(MemoryPool)];
static MemoryPool *defaultMemPool = (MemoryPool *)defaultMemPool_space;
const size_t MemoryPool::defaultGranularity;
// zero-initialized
MallocMutex  MemoryPool::memPoolListLock;
// TODO: move huge page status to default pool, because that's its states
HugePagesStatus hugePages;

// Slab block is 16KB-aligned. To prvent false sharing, separate locally-accessed
// fields and fields commonly accessed by not owner threads.
class GlobalBlockFields : public BlockI {
protected:
    FreeObject  *publicFreeList;
    Block       *nextPrivatizable;
};

class LocalBlockFields : public GlobalBlockFields {
protected:
    size_t       __pad_local_fields[(blockHeaderAlignment -
                                     sizeof(GlobalBlockFields))/sizeof(size_t)];

    Block       *next;
    Block       *previous;        /* Use double linked list to speed up removal */
    uint16_t     objectSize;
    ThreadId     owner;
    FreeObject  *bumpPtr;         /* Bump pointer moves from the end to the beginning of a block */
    FreeObject  *freeList;
    BackRefIdx   backRefIdx;
    unsigned int allocatedCount;  /* Number of objects allocated (obviously by the owning thread) */
    bool         isFull;
    bool         orphaned;

    friend void *BootStrapBlocks::allocate(MemoryPool *memPool, size_t size);
    friend class FreeBlockPool;
    friend class StartupBlock;
    friend class LifoList;
    friend Block *MemoryPool::getEmptyBlock(size_t size);
};

class Block : public LocalBlockFields {
    size_t       __pad_public_fields[(2*blockHeaderAlignment -
                                      sizeof(LocalBlockFields))/sizeof(size_t)];
public:
    bool empty() const { return allocatedCount==0 && publicFreeList==NULL; }
    inline FreeObject* allocate();
    inline FreeObject *allocateFromFreeList();
    inline bool emptyEnoughToUse();
    bool freeListNonNull() { return freeList; }
    void freePublicObject(FreeObject *objectToFree);
    inline void freeOwnObject(MemoryPool *memPool, TLSData *tls, void *object);
    void makeEmpty();
    void privatizePublicFreeList();
    void restoreBumpPtr();
    void privatizeOrphaned(Bin *bin);
    void shareOrphaned(const Bin *bin);
    unsigned int getSize() const {
        MALLOC_ASSERT(isStartupAllocObject() || objectSize<minLargeObjectSize,
                      "Invalid object size");
        return objectSize;
    }
    const BackRefIdx *getBackRefIdx() const { return &backRefIdx; }
    bool ownBlock() const { return !orphaned && owner.own(); }
    bool isStartupAllocObject() const { return objectSize == startupAllocObjSizeMark; }
    inline FreeObject *findObjectToFree(void *object) const;
    bool checkFreePrecond(void *object) const {
        if (allocatedCount>0) {
            if (startupAllocObjSizeMark == objectSize) // startup block
                return object<=bumpPtr;
            else
                return allocatedCount <= (slabSize-sizeof(Block))/objectSize
                       && (!bumpPtr || object>bumpPtr);
        }
        return false;
    }
    const BackRefIdx *getBackRef() const { return &backRefIdx; }
    void initEmptyBlock(Bin* tlsBin, size_t size);

protected:
    void cleanBlockHeader();

private:
    static const float emptyEnoughRatio; /* "Reactivate" a block if this share of its objects is free. */

    inline FreeObject *allocateFromBumpPtr();
    inline FreeObject *findAllocatedObject(const void *address) const;
    inline bool isProperlyPlaced(const void *object) const;

    friend class Bin;
    friend class TLSData;
    friend void MemoryPool::destroy();
};

const float Block::emptyEnoughRatio = 1.0 / 4.0;

class Bin {
    Block      *activeBlk;
    Block      *mailbox;
    MallocMutex mailLock;

public:
    inline Block* getActiveBlock() const { return activeBlk; }
    void resetActiveBlock() { activeBlk = 0; }
    bool activeBlockUnused() const { return activeBlk && !activeBlk->allocatedCount; }
    inline void setActiveBlock(Block *block);
    inline Block* setPreviousBlockActive();
    Block* getPublicFreeListBlock();
    void moveBlockToBinFront(Block *block);
    void processLessUsedBlock(MemoryPool *memPool, Block *block);

    void outofTLSBin (Block* block);
    void verifyTLSBin (size_t size) const;
    void pushTLSBin(Block* block);

    void verifyInitState() const {
        MALLOC_ASSERT( activeBlk == 0, ASSERT_TEXT );
        MALLOC_ASSERT( mailbox == 0, ASSERT_TEXT );
    }

    friend void Block::freePublicObject (FreeObject *objectToFree);
};

/********* End of the data structures                    **************/

/*
 * There are bins for all 8 byte aligned objects less than this segregated size; 8 bins in total
 */
const uint32_t minSmallObjectIndex = 0;
const uint32_t numSmallObjectBins = 8;
const uint32_t maxSmallObjectSize = 64;

/*
 * There are 4 bins between each couple of powers of 2 [64-128-256-...]
 * from maxSmallObjectSize till this size; 16 bins in total
 */
const uint32_t minSegregatedObjectIndex = minSmallObjectIndex+numSmallObjectBins;
const uint32_t numSegregatedObjectBins = 16;
const uint32_t maxSegregatedObjectSize = 1024;

/*
 * And there are 5 bins with allocation sizes that are multiples of estimatedCacheLineSize
 * and selected to fit 9, 6, 4, 3, and 2 allocations in a block.
 */
const uint32_t minFittingIndex = minSegregatedObjectIndex+numSegregatedObjectBins;
const uint32_t numFittingBins = 5;

const uint32_t fittingAlignment = estimatedCacheLineSize;

#define SET_FITTING_SIZE(N) ( (slabSize-sizeof(Block))/N ) & ~(fittingAlignment-1)
// For blockSize=16*1024, sizeof(Block)=2*estimatedCacheLineSize and fittingAlignment=estimatedCacheLineSize,
// the comments show the fitting sizes and the amounts left unused for estimatedCacheLineSize=64/128:
const uint32_t fittingSize1 = SET_FITTING_SIZE(9); // 1792/1792 128/000
const uint32_t fittingSize2 = SET_FITTING_SIZE(6); // 2688/2688 128/000
const uint32_t fittingSize3 = SET_FITTING_SIZE(4); // 4032/3968 128/256
const uint32_t fittingSize4 = SET_FITTING_SIZE(3); // 5376/5376 128/000
const uint32_t fittingSize5 = SET_FITTING_SIZE(2); // 8128/8064 000/000
#undef SET_FITTING_SIZE

/*
 * The total number of thread-specific Block-based bins
 */
const uint32_t numBlockBins = minFittingIndex+numFittingBins;

/*
 * Objects of this size and larger are considered large objects.
 */
const uint32_t minLargeObjectSize = fittingSize5 + 1;

/*
 * Default granularity of memory pools
 */

#if USE_WINTHREAD
const size_t scalableMallocPoolGranularity = 64*1024; // for VirtualAlloc use
#else
const size_t scalableMallocPoolGranularity = 4*1024;  // page size, for mmap use
#endif

/*
 * Per-thread pool of slab blocks. Idea behind it is to not share with other
 * threads memory that are likely in local cache(s) of our CPU.
 */
class FreeBlockPool {
    Block      *head;
    Block      *tail;
    int         size;
    Backend    *backend;
    bool        lastAccessMiss;
    void insertBlock(Block *block);
public:
    static const int POOL_HIGH_MARK = 32;
    static const int POOL_LOW_MARK  = 8;

    class ResOfGet {
        ResOfGet();
    public:
        Block* block;
        bool   lastAccMiss;
        ResOfGet(Block *b, bool lastMiss) : block(b), lastAccMiss(lastMiss) {}
    };

    // allocated in zero-initialized memory
    FreeBlockPool(Backend *bknd) : backend(bknd) {}
    ResOfGet getBlock();
    void returnBlock(Block *block);
    bool releaseAllBlocks();
};

template<int LOW_MARK, int HIGH_MARK>
class LocalLOC {
    static const size_t MAX_TOTAL_SIZE = 4*1024*1024;

    // TODO: can single-linked list be faster here?
    LargeMemoryBlock *head,
                     *tail;
    intptr_t          lastSeenOSCallsCnt,
                      lastUsedOSCallsCnt;
    size_t            totalSize;
    int               numOfBlocks;
public:
    bool put(LargeMemoryBlock *object, ExtMemoryPool *extMemPool);
    LargeMemoryBlock *get(size_t size);
    bool clean(ExtMemoryPool *extMemPool);
    void allocatorCalledHook(ExtMemoryPool *extMemPool);
#if __TBB_MALLOC_WHITEBOX_TEST
    LocalLOC() : head(NULL), tail(NULL), lastSeenOSCallsCnt(0),
                 lastUsedOSCallsCnt(0), totalSize(0),
                 numOfBlocks(0) {}
    static size_t getMaxSize() { return MAX_TOTAL_SIZE; }
#else
    // no ctor, object must be created in zero-initialized memory
#endif
};

class TLSData {
#if USE_PTHREAD
    MemoryPool   *memPool;
#endif
public:
    Bin           bin[numBlockBinLimit];
    FreeBlockPool freeSlabBlocks;
    LocalLOC<8,32> lloc;
#if USE_PTHREAD
    TLSData(MemoryPool *mPool, Backend *bknd) : memPool(mPool), freeSlabBlocks(bknd) {}
    MemoryPool *getMemPool() const { return memPool; }
#else
    TLSData(MemoryPool * /*memPool*/, Backend *bknd) : freeSlabBlocks(bknd) {}
#endif
    void release(MemoryPool *mPool);
};

TLSData *TLSKey::createTLS(MemoryPool *memPool, Backend *backend)
{
    MALLOC_ASSERT( sizeof(TLSData) >= sizeof(Bin) * numBlockBins + sizeof(FreeBlockPool), ASSERT_TEXT );
    TLSData* tls = (TLSData*) memPool->bootStrapBlocks.allocate(memPool, sizeof(TLSData));
    if ( !tls )
        return NULL;
    new(tls) TLSData(memPool, backend);
    /* the block contains zeroes after bootStrapMalloc, so bins are initialized */
#if MALLOC_DEBUG
    for (uint32_t i = 0; i < numBlockBinLimit; i++)
        tls->bin[i].verifyInitState();
#endif
    setThreadMallocTLS(tls);
    return tls;
}

bool ExtMemoryPool::releaseTLCaches()
{
    bool released = false;

    if (TLSData *tlsData = tlsPointerKey.getThreadMallocTLS()) {
        released = tlsData->freeSlabBlocks.releaseAllBlocks();
        released |= tlsData->lloc.clean(this);

        // active blocks can be not used, so return them to backend
        for (uint32_t i=0; i<numBlockBinLimit; i++)
            if (tlsData->bin[i].activeBlockUnused()) {
                Block *block = tlsData->bin[i].getActiveBlock();
                tlsData->bin[i].outofTLSBin(block);
                // slab blocks in user's pools do not have valid backRefIdx
                if (!userPool())
                    removeBackRef(*(block->getBackRefIdx()));
                backend.putSlabBlock(block);

                released = true;
            }
    }
    return released;
}


#if MALLOC_CHECK_RECURSION
MallocMutex RecursiveMallocCallProtector::rmc_mutex;
pthread_t   RecursiveMallocCallProtector::owner_thread;
void       *RecursiveMallocCallProtector::autoObjPtr;
bool        RecursiveMallocCallProtector::mallocRecursionDetected;
#if __FreeBSD__
bool        RecursiveMallocCallProtector::canUsePthread;
#endif

#endif

/*********** End code to provide thread ID and a TLS pointer **********/

static void *internalMalloc(size_t size);
static void internalFree(void *object);
static void *internalPoolMalloc(MemoryPool* mPool, size_t size);
static bool internalPoolFree(MemoryPool *mPool, void *object);

#if !MALLOC_DEBUG
#if __INTEL_COMPILER || _MSC_VER
#define NOINLINE(decl) __declspec(noinline) decl
#define ALWAYSINLINE(decl) __forceinline decl
#elif __GNUC__
#define NOINLINE(decl) decl __attribute__ ((noinline))
#define ALWAYSINLINE(decl) decl __attribute__ ((always_inline))
#else
#define NOINLINE(decl) decl
#define ALWAYSINLINE(decl) decl
#endif

static NOINLINE( void doInitialization() );
ALWAYSINLINE( bool isMallocInitialized() );

#undef ALWAYSINLINE
#undef NOINLINE
#endif /* !MALLOC_DEBUG */


/********* Now some rough utility code to deal with indexing the size bins. **************/

/*
 * Given a number return the highest non-zero bit in it. It is intended to work with 32-bit values only.
 * Moreover, on IPF, for sake of simplicity and performance, it is narrowed to only serve for 64 to 1023.
 * This is enough for current algorithm of distribution of sizes among bins.
 * __TBB_Log2 is not used here to minimize dependencies on TBB specific sources.
 */
#if _WIN64 && _MSC_VER>=1400 && !__INTEL_COMPILER
extern "C" unsigned char _BitScanReverse( unsigned long* i, unsigned long w );
#pragma intrinsic(_BitScanReverse)
#endif
static inline unsigned int highestBitPos(unsigned int n)
{
    MALLOC_ASSERT( n>=64 && n<1024, ASSERT_TEXT ); // only needed for bsr array lookup, but always true
    unsigned int pos;
#if __ARCH_x86_32||__ARCH_x86_64

# if __linux__||__APPLE__||__FreeBSD__||__NetBSD__||__sun||__MINGW32__
    __asm__ ("bsr %1,%0" : "=r"(pos) : "r"(n));
# elif (_WIN32 && (!_WIN64 || __INTEL_COMPILER))
    __asm
    {
        bsr eax, n
        mov pos, eax
    }
# elif _WIN64 && _MSC_VER>=1400
    _BitScanReverse((unsigned long*)&pos, (unsigned long)n);
# else
#   error highestBitPos() not implemented for this platform
# endif
#elif __arm__
    __asm__ __volatile__
    (
       "clz %0, %1\n"
       "rsb %0, %0, %2\n"
       :"=r" (pos) :"r" (n), "I" (31)
    );
#else
    static unsigned int bsr[16] = {0/*N/A*/,6,7,7,8,8,8,8,9,9,9,9,9,9,9,9};
    pos = bsr[ n>>6 ];
#endif /* __ARCH_* */
    return pos;
}

/*
 * Depending on indexRequest, for a given size return either the index into the bin
 * for objects of this size, or the actual size of objects in this bin.
 */
template<bool indexRequest>
static unsigned int getIndexOrObjectSize (unsigned int size)
{
    if (size <= maxSmallObjectSize) { // selection from 4/8/16/24/32/40/48/56/64
         /* Index 0 holds up to 8 bytes, Index 1 16 and so forth */
        return indexRequest ? (size - 1) >> 3 : alignUp(size,8);
    }
    else if (size <= maxSegregatedObjectSize ) { // 80/96/112/128 / 160/192/224/256 / 320/384/448/512 / 640/768/896/1024
        unsigned int order = highestBitPos(size-1); // which group of bin sizes?
        MALLOC_ASSERT( 6<=order && order<=9, ASSERT_TEXT );
        if (indexRequest)
            return minSegregatedObjectIndex - (4*6) - 4 + (4*order) + ((size-1)>>(order-2));
        else {
            unsigned int alignment = 128 >> (9-order); // alignment in the group
            MALLOC_ASSERT( alignment==16 || alignment==32 || alignment==64 || alignment==128, ASSERT_TEXT );
            return alignUp(size,alignment);
        }
    }
    else {
        if( size <= fittingSize3 ) {
            if( size <= fittingSize2 ) {
                if( size <= fittingSize1 )
                    return indexRequest ? minFittingIndex : fittingSize1;
                else
                    return indexRequest ? minFittingIndex+1 : fittingSize2;
            } else
                return indexRequest ? minFittingIndex+2 : fittingSize3;
        } else {
            if( size <= fittingSize5 ) {
                if( size <= fittingSize4 )
                    return indexRequest ? minFittingIndex+3 : fittingSize4;
                else
                    return indexRequest ? minFittingIndex+4 : fittingSize5;
            } else {
                MALLOC_ASSERT( 0,ASSERT_TEXT ); // this should not happen
                return ~0U;
            }
        }
    }
}

static unsigned int getIndex (unsigned int size)
{
    return getIndexOrObjectSize</*indexRequest*/true>(size);
}

static unsigned int getObjectSize (unsigned int size)
{
    return getIndexOrObjectSize</*indexRequest*/false>(size);
}


void *BootStrapBlocks::allocate(MemoryPool *memPool, size_t size)
{
    FreeObject *result;

    MALLOC_ASSERT( size == sizeof(TLSData), ASSERT_TEXT );

    { // Lock with acquire
        MallocMutex::scoped_lock scoped_cs(bootStrapLock);

        if( bootStrapObjectList) {
            result = bootStrapObjectList;
            bootStrapObjectList = bootStrapObjectList->next;
        } else {
            if (!bootStrapBlock) {
                bootStrapBlock = memPool->getEmptyBlock(size);
                if (!bootStrapBlock) return NULL;
            }
            result = bootStrapBlock->bumpPtr;
            bootStrapBlock->bumpPtr = (FreeObject *)((uintptr_t)bootStrapBlock->bumpPtr - bootStrapBlock->objectSize);
            if ((uintptr_t)bootStrapBlock->bumpPtr < (uintptr_t)bootStrapBlock+sizeof(Block)) {
                bootStrapBlock->bumpPtr = NULL;
                bootStrapBlock->next = bootStrapBlockUsed;
                bootStrapBlockUsed = bootStrapBlock;
                bootStrapBlock = NULL;
            }
        }
    } // Unlock with release

    memset (result, 0, size);
    return (void*)result;
}

void BootStrapBlocks::free(void* ptr)
{
    MALLOC_ASSERT( ptr, ASSERT_TEXT );
    { // Lock with acquire
        MallocMutex::scoped_lock scoped_cs(bootStrapLock);
        ((FreeObject*)ptr)->next = bootStrapObjectList;
        bootStrapObjectList = (FreeObject*)ptr;
    } // Unlock with release
}

void BootStrapBlocks::reset()
{
    bootStrapBlock = bootStrapBlockUsed = NULL;
    bootStrapObjectList = NULL;
}

#if !(FREELIST_NONBLOCKING)
static MallocMutex publicFreeListLock; // lock for changes of publicFreeList
#endif

const uintptr_t UNUSABLE = 0x1;
inline bool isSolidPtr( void* ptr )
{
    return (UNUSABLE|(uintptr_t)ptr)!=UNUSABLE;
}
inline bool isNotForUse( void* ptr )
{
    return (uintptr_t)ptr==UNUSABLE;
}

/********* End rough utility code  **************/

#ifdef FINE_GRAIN_LOCKS
/* LifoList assumes zero initialization so a vector of it can be created
 * by just allocating some space with no call to constructor.
 * On Linux, it seems to be necessary to avoid linking with C++ libraries.
 *
 * By usage convention there is no race on the initialization. */
LifoList::LifoList( ) : top(NULL)
{
    // MallocMutex assumes zero initialization
    memset(&lock, 0, sizeof(MallocMutex));
}

void LifoList::push(Block *block)
{
    MallocMutex::scoped_lock scoped_cs(lock);
    block->next = top;
    top = block;
}

Block *LifoList::pop()
{
    Block *block=NULL;
    if (!top) goto done;
    {
        MallocMutex::scoped_lock scoped_cs(lock);
        if (!top) goto done;
        block = top;
        top = block->next;
    }
done:
    return block;
}

#endif /* FINE_GRAIN_LOCKS     */

/********* Thread and block related code      *************/

TLSData* MemoryPool::getTLS(bool create)
{
    TLSData* tls = extMemPool.tlsPointerKey.getThreadMallocTLS();
    if( create && !tls ) {
        tls = extMemPool.tlsPointerKey.createTLS(this, &extMemPool.backend);
        MALLOC_ASSERT( tls, ASSERT_TEXT );
    }
    return tls;
}

/*
 * Return the bin for the given size.
 */
Bin* MemoryPool::getAllocationBin(TLSData* tls, size_t size)
{
    return tls->bin + getIndex(size);
}

/* Return an empty uninitialized block in a non-blocking fashion. */
Block *MemoryPool::getEmptyBlock(size_t size)
{
    FreeBlockPool::ResOfGet resOfGet(NULL, false);
    Block *result = NULL, *b;
    TLSData* tls = extMemPool.tlsPointerKey.getThreadMallocTLS();

    if (tls)
        resOfGet = tls->freeSlabBlocks.getBlock();
    if (resOfGet.block) {
        result = resOfGet.block;
    } else {
        int i, num = resOfGet.lastAccMiss? Backend::numOfSlabAllocOnMiss : 1;
        BackRefIdx backRefIdx[Backend::numOfSlabAllocOnMiss];

        result = static_cast<Block*>(extMemPool.backend.getSlabBlock(num));
        if (!result) return NULL;

        if (!extMemPool.userPool())
            for (i=0; i<num; i++) {
                backRefIdx[i] = BackRefIdx::newBackRef(/*largeObj=*/false);
                if (backRefIdx[i].isInvalid()) {
                    // roll back resource allocation
                    for (int j=0; j<i; j++)
                        removeBackRef(backRefIdx[j]);
                    Block *b;
                    for (b=result, i=0; i<num;
                         b=(Block*)((uintptr_t)b+slabSize), i++)
                        extMemPool.backend.putSlabBlock(b);
                    return NULL;
                }
            }
        // resources were allocated, register blocks
        for (b=result, i=0; i<num; b=(Block*)((uintptr_t)b+slabSize), i++) {
            // slab block in user's pool must have invalid backRefIdx
            if (extMemPool.userPool()) {
                new (&b->backRefIdx) BackRefIdx();
            } else {
                setBackRef(backRefIdx[i], b);
                b->backRefIdx = backRefIdx[i];
            }
            // all but first one go to per-thread pool
            if (i > 0) {
                MALLOC_ASSERT(tls, ASSERT_TEXT);
                tls->freeSlabBlocks.returnBlock(b);
            }
        }
    }
    if (result) {
        result->initEmptyBlock(tls? tls->bin : NULL, size);
        STAT_increment(result->owner, getIndex(result->objectSize), allocBlockNew);
    }
    return result;
}

void MemoryPool::returnEmptyBlock(Block *block, bool poolTheBlock)
{
    block->makeEmpty();
    if (poolTheBlock) {
        extMemPool.tlsPointerKey.getThreadMallocTLS()->freeSlabBlocks.returnBlock(block);
    }
    else {
        // slab blocks in user's pools do not have valid backRefIdx
        if (!extMemPool.userPool())
            removeBackRef(*(block->getBackRefIdx()));
        extMemPool.backend.putSlabBlock(block);
    }
}

bool ExtMemoryPool::init(intptr_t poolId, rawAllocType rawAlloc,
                         rawFreeType rawFree, size_t granularity,
                         bool keepAllMemory, bool fixedPool)
{
    this->poolId = poolId;
    this->rawAlloc = rawAlloc;
    this->rawFree = rawFree;
    this->granularity = granularity;
    this->keepAllMemory = keepAllMemory;
    this->fixedPool = fixedPool;
    this->delayRegsReleasing = false;
    initTLS();
    // allocate initial region for user's objects placement
    return backend.bootstrap(this);
}

void ExtMemoryPool::initTLS() { new (&tlsPointerKey) TLSKey(); }

bool MemoryPool::init(intptr_t poolId, const MemPoolPolicy *policy)
{
    if (!extMemPool.init(poolId, policy->pAlloc, policy->pFree,
               policy->granularity? policy->granularity : defaultGranularity,
               policy->keepAllMemory, policy->fixedPool))
        return false;
    {
        MallocMutex::scoped_lock lock(memPoolListLock);
        next = defaultMemPool->next;
        defaultMemPool->next = this;
        prev = defaultMemPool;
        if (next)
            next->prev = this;
    }
    return true;
}

void MemoryPool::reset()
{
    // memory is not releasing during pool reset
    // TODO: mark regions to release unused on next reset()
    extMemPool.delayRegionsReleasing(true);

    bootStrapBlocks.reset();
    orphanedBlocks.reset();
    extMemPool.reset();

    extMemPool.initTLS();
    extMemPool.delayRegionsReleasing(false);
}

void MemoryPool::destroy()
{
    {
        MallocMutex::scoped_lock lock(memPoolListLock);
        // remove itself from global pool list
        if (prev)
            prev->next = next;
        if (next)
            next->prev = prev;
    }
    // slab blocks in non-default pool do not have backreferencies,
    // only large objects do
    for (LargeMemoryBlock *lmb = extMemPool.lmbList.getHead(); lmb; ) {
        LargeMemoryBlock *next = lmb->gNext;
        if (extMemPool.userPool())
            removeBackRef(lmb->backRefIdx);
        lmb = next;
    }
    extMemPool.destroy();
}

void MemoryPool::processThreadShutdown(TLSData *tlsData)
{
    tlsData->release(this);
    bootStrapBlocks.free(tlsData);
    clearTLS();
}

void Bin::verifyTLSBin (size_t size) const
{
    suppress_unused_warning(size);
#if MALLOC_DEBUG
/* The debug version verifies the TLSBin as needed */
    uint32_t objSize = getObjectSize(size);

    if (activeBlk) {
        MALLOC_ASSERT( activeBlk->owner.own(), ASSERT_TEXT );
        MALLOC_ASSERT( activeBlk->objectSize == objSize, ASSERT_TEXT );
#if MALLOC_DEBUG>1
        for (Block* temp = activeBlk->next; temp; temp=temp->next) {
            MALLOC_ASSERT( temp!=activeBlk, ASSERT_TEXT );
            MALLOC_ASSERT( temp->owner.own(), ASSERT_TEXT );
            MALLOC_ASSERT( temp->objectSize == objSize, ASSERT_TEXT );
            MALLOC_ASSERT( temp->previous->next == temp, ASSERT_TEXT );
            if (temp->next) {
                MALLOC_ASSERT( temp->next->previous == temp, ASSERT_TEXT );
            }
        }
        for (Block* temp = activeBlk->previous; temp; temp=temp->previous) {
            MALLOC_ASSERT( temp!=activeBlk, ASSERT_TEXT );
            MALLOC_ASSERT( temp->owner.own(), ASSERT_TEXT );
            MALLOC_ASSERT( temp->objectSize == objSize, ASSERT_TEXT );
            MALLOC_ASSERT( temp->next->previous == temp, ASSERT_TEXT );
            if (temp->previous) {
                MALLOC_ASSERT( temp->previous->next == temp, ASSERT_TEXT );
            }
        }
#endif /* MALLOC_DEBUG>1 */
    }
#endif /* MALLOC_DEBUG */
}

/*
 * Add a block to the start of this tls bin list.
 */
void Bin::pushTLSBin(Block* block)
{
    /* The objectSize should be defined and not a parameter
       because the function is applied to partially filled blocks as well */
    unsigned int size = block->objectSize;

    MALLOC_ASSERT( block->owner == ThreadId::get(), ASSERT_TEXT );
    MALLOC_ASSERT( block->objectSize != 0, ASSERT_TEXT );
    MALLOC_ASSERT( block->next == NULL, ASSERT_TEXT );
    MALLOC_ASSERT( block->previous == NULL, ASSERT_TEXT );

    MALLOC_ASSERT( this, ASSERT_TEXT );
    verifyTLSBin(size);

    block->next = activeBlk;
    if( activeBlk ) {
        block->previous = activeBlk->previous;
        activeBlk->previous = block;
        if( block->previous )
            block->previous->next = block;
    } else {
        activeBlk = block;
    }

    verifyTLSBin(size);
}

/*
 * Take a block out of its tls bin (e.g. before removal).
 */
void Bin::outofTLSBin(Block* block)
{
    unsigned int size = block->objectSize;

    MALLOC_ASSERT( block->owner == ThreadId::get(), ASSERT_TEXT );
    MALLOC_ASSERT( block->objectSize != 0, ASSERT_TEXT );

    MALLOC_ASSERT( this, ASSERT_TEXT );
    verifyTLSBin(size);

    if (block == activeBlk) {
        activeBlk = block->previous? block->previous : block->next;
    }
    /* Delink the block */
    if (block->previous) {
        MALLOC_ASSERT( block->previous->next == block, ASSERT_TEXT );
        block->previous->next = block->next;
    }
    if (block->next) {
        MALLOC_ASSERT( block->next->previous == block, ASSERT_TEXT );
        block->next->previous = block->previous;
    }
    block->next = NULL;
    block->previous = NULL;

    verifyTLSBin(size);
}

Block* Bin::getPublicFreeListBlock()
{
    Block* block;
    MALLOC_ASSERT( this, ASSERT_TEXT );
    // if this method is called, active block usage must be unsuccesful
    MALLOC_ASSERT( !activeBlk && !mailbox || activeBlk && activeBlk->isFull, ASSERT_TEXT );

// the counter should be changed    STAT_increment(getThreadId(), ThreadCommonCounters, lockPublicFreeList);
    {
        MallocMutex::scoped_lock scoped_cs(mailLock);
        block = mailbox;
        if( block ) {
            MALLOC_ASSERT( block->ownBlock(), ASSERT_TEXT );
            MALLOC_ASSERT( !isNotForUse(block->nextPrivatizable), ASSERT_TEXT );
            mailbox = block->nextPrivatizable;
            block->nextPrivatizable = (Block*) this;
        }
    }
    if( block ) {
        MALLOC_ASSERT( isSolidPtr(block->publicFreeList), ASSERT_TEXT );
        block->privatizePublicFreeList();
    }
    return block;
}

bool Block::emptyEnoughToUse()
{
    const float threshold = (slabSize - sizeof(Block)) * (1-emptyEnoughRatio);

    if (bumpPtr) {
        /* If we are still using a bump ptr for this block it is empty enough to use. */
        STAT_increment(owner, getIndex(objectSize), examineEmptyEnough);
        isFull = false;
        return 1;
    }

    /* allocatedCount shows how many objects in the block are in use; however it still counts
       blocks freed by other threads; so prior call to privatizePublicFreeList() is recommended */
    isFull = (allocatedCount*objectSize > threshold)? true: false;
#if COLLECT_STATISTICS
    if (isFull)
        STAT_increment(owner, getIndex(objectSize), examineNotEmpty);
    else
        STAT_increment(owner, getIndex(objectSize), examineEmptyEnough);
#endif
    return !isFull;
}

/* Restore the bump pointer for an empty block that is planned to use */
void Block::restoreBumpPtr()
{
    MALLOC_ASSERT( allocatedCount == 0, ASSERT_TEXT );
    MALLOC_ASSERT( publicFreeList == NULL, ASSERT_TEXT );
    STAT_increment(owner, getIndex(objectSize), freeRestoreBumpPtr);
    bumpPtr = (FreeObject *)((uintptr_t)this + slabSize - objectSize);
    freeList = NULL;
    isFull = 0;
}

void Block::freeOwnObject(MemoryPool *memPool, TLSData *tls, void *object)
{
    allocatedCount--;
    MALLOC_ASSERT( allocatedCount < (slabSize-sizeof(Block))/objectSize, ASSERT_TEXT );
#if COLLECT_STATISTICS
    if (getActiveBlock(memPool->getAllocationBin(block->objectSize)) != block)
        STAT_increment(myTid, getIndex(block->objectSize), freeToInactiveBlock);
    else
        STAT_increment(myTid, getIndex(block->objectSize), freeToActiveBlock);
#endif
    if (allocatedCount==0 && publicFreeList==NULL) {
        // The bump pointer is about to be restored for the block,
        // no need to find objectToFree here (this is costly).

        // if the last object of a slab is freed, the slab cannot be marked full
        MALLOC_ASSERT(!isFull, ASSERT_TEXT);
        memPool->getAllocationBin(tls, objectSize)->
            processLessUsedBlock(memPool, this);
    } else {
        FreeObject *objectToFree = findObjectToFree(object);
        objectToFree->next = freeList;
        freeList = objectToFree;

        if (isFull) {
            if (emptyEnoughToUse())
                memPool->getAllocationBin(tls, objectSize)->moveBlockToBinFront(this);
        }
    }
}

void Block::freePublicObject (FreeObject *objectToFree)
{
    FreeObject *localPublicFreeList;

    MALLOC_ITT_SYNC_RELEASING(&publicFreeList);
#if FREELIST_NONBLOCKING
    FreeObject *temp = publicFreeList;
    do {
        localPublicFreeList = objectToFree->next = temp;
        temp = (FreeObject*)AtomicCompareExchange(
                                (intptr_t&)publicFreeList,
                                (intptr_t)objectToFree, (intptr_t)localPublicFreeList );
        // no backoff necessary because trying to make change, not waiting for a change
    } while( temp != localPublicFreeList );
#else
    STAT_increment(getThreadId(), ThreadCommonCounters, lockPublicFreeList);
    {
        MallocMutex::scoped_lock scoped_cs(publicFreeListLock);
        localPublicFreeList = objectToFree->next = publicFreeList;
        publicFreeList = objectToFree;
    }
#endif

    if( localPublicFreeList==NULL ) {
        // if the block is abandoned, its nextPrivatizable pointer should be UNUSABLE
        // otherwise, it should point to the bin the block belongs to.
        // reading nextPrivatizable is thread-safe below, because:
        // 1) the executing thread atomically got localPublicFreeList==NULL and changed it to non-NULL;
        // 2) only owning thread can change it back to NULL,
        // 3) but it can not be done until the block is put to the mailbox
        // So the executing thread is now the only one that can change nextPrivatizable
        if( !isNotForUse(nextPrivatizable) ) {
            MALLOC_ASSERT( nextPrivatizable!=NULL, ASSERT_TEXT );
            MALLOC_ASSERT( owner.defined(), ASSERT_TEXT );
            Bin* theBin = (Bin*) nextPrivatizable;
            MallocMutex::scoped_lock scoped_cs(theBin->mailLock);
            nextPrivatizable = theBin->mailbox;
            theBin->mailbox = this;
        } else {
            MALLOC_ASSERT( !owner.defined(), ASSERT_TEXT );
        }
    }
    STAT_increment(ThreadId::get(), ThreadCommonCounters, freeToOtherThread);
    STAT_increment(owner, getIndex(objectSize), freeByOtherThread);
}

void Block::privatizePublicFreeList()
{
    FreeObject *temp, *localPublicFreeList;

    MALLOC_ASSERT( owner.own(), ASSERT_TEXT );
#if FREELIST_NONBLOCKING
    temp = publicFreeList;
    do {
        localPublicFreeList = temp;
        temp = (FreeObject*)AtomicCompareExchange(
                                (intptr_t&)publicFreeList,
                                0, (intptr_t)localPublicFreeList);
        // no backoff necessary because trying to make change, not waiting for a change
    } while( temp != localPublicFreeList );
#else
    STAT_increment(owner, ThreadCommonCounters, lockPublicFreeList);
    {
        MallocMutex::scoped_lock scoped_cs(publicFreeListLock);
        localPublicFreeList = publicFreeList;
        publicFreeList = NULL;
    }
    temp = localPublicFreeList;
#endif
    MALLOC_ITT_SYNC_ACQUIRED(&publicFreeList);

    MALLOC_ASSERT( localPublicFreeList && localPublicFreeList==temp, ASSERT_TEXT ); // there should be something in publicFreeList!
    if( !isNotForUse(temp) ) { // return/getPartialBlock could set it to UNUSABLE
        MALLOC_ASSERT( allocatedCount <= (slabSize-sizeof(Block))/objectSize, ASSERT_TEXT );
        /* other threads did not change the counter freeing our blocks */
        allocatedCount--;
        while( isSolidPtr(temp->next) ){ // the list will end with either NULL or UNUSABLE
            temp = temp->next;
            allocatedCount--;
        }
        MALLOC_ASSERT( allocatedCount < (slabSize-sizeof(Block))/objectSize, ASSERT_TEXT );
        /* merge with local freeList */
        temp->next = freeList;
        freeList = localPublicFreeList;
        STAT_increment(owner, getIndex(objectSize), allocPrivatized);
    }
}

void Block::privatizeOrphaned(Bin* bin)
{
    next = NULL;
    previous = NULL;
    MALLOC_ASSERT( publicFreeList!=NULL, ASSERT_TEXT );
    /* There is not a race here since no other thread owns this block */
    MALLOC_ASSERT( !owner.defined(), ASSERT_TEXT );
    owner = ThreadId::get();
    MALLOC_ASSERT(orphaned, ASSERT_TEXT);
    orphaned = false;
    // It is safe to change nextPrivatizable, as publicFreeList is not null
    MALLOC_ASSERT( isNotForUse(nextPrivatizable), ASSERT_TEXT );
    nextPrivatizable = (Block*)bin;
    // the next call is required to change publicFreeList to 0
    privatizePublicFreeList();
    if( allocatedCount ) {
        emptyEnoughToUse(); // check its fullness and set result->isFull
    } else {
        restoreBumpPtr();
    }
    MALLOC_ASSERT( !isNotForUse(publicFreeList), ASSERT_TEXT );
}

void Block::shareOrphaned(const Bin *bin)
{
    MALLOC_ASSERT( bin, ASSERT_TEXT );
    STAT_increment(owner, index, freeBlockPublic);
    MALLOC_ASSERT(!orphaned, ASSERT_TEXT);
    orphaned = true;
    // need to set publicFreeList to non-zero, so other threads
    // will not change nextPrivatizable and it can be zeroed.
    if ((intptr_t)nextPrivatizable==(intptr_t)bin) {
        void* oldval;
#if FREELIST_NONBLOCKING
        oldval = (void*)AtomicCompareExchange((intptr_t&)publicFreeList, (intptr_t)UNUSABLE, 0);
#else
        STAT_increment(owner, ThreadCommonCounters, lockPublicFreeList);
        {
            MallocMutex::scoped_lock scoped_cs(publicFreeListLock);
            if ( (oldval=publicFreeList)==NULL )
                (uintptr_t&)(publicFreeList) = UNUSABLE;
        }
#endif
        if ( oldval!=NULL ) {
            // another thread freed an object; we need to wait until it finishes.
            // I believe there is no need for exponential backoff, as the wait here is not for a lock;
            // but need to yield, so the thread we wait has a chance to run.
            int count = 256;
            while( (intptr_t)const_cast<Block* volatile &>(nextPrivatizable)==(intptr_t)bin ) {
                if (--count==0) {
                    do_yield();
                    count = 256;
                }
            }
        }
    } else {
        MALLOC_ASSERT( isSolidPtr(publicFreeList), ASSERT_TEXT );
    }
    MALLOC_ASSERT( publicFreeList!=NULL, ASSERT_TEXT );
    // now it is safe to change our data
    previous = NULL;
    owner.undef();
    // it is caller responsibility to ensure that the list of blocks
    // formed by nextPrivatizable pointers is kept consistent if required.
    // if only called from thread shutdown code, it does not matter.
    (uintptr_t&)(nextPrivatizable) = UNUSABLE;
}

void Block::cleanBlockHeader()
{
    next = NULL;
    previous = NULL;
    freeList = NULL;
    allocatedCount = 0;
    isFull = 0;
    orphaned = false;

    publicFreeList = NULL;
}

void Block::initEmptyBlock(Bin* tlsBin, size_t size)
{
    // Having getIndex and getObjectSize called next to each other
    // allows better compiler optimization as they basically share the code.
    unsigned int index = getIndex(size);
    unsigned int objSz = getObjectSize(size);

    cleanBlockHeader();
    objectSize = objSz;
    owner = ThreadId::get();
    // bump pointer should be prepared for first allocation - thus mode it down to objectSize
    bumpPtr = (FreeObject *)((uintptr_t)this + slabSize - objectSize);

    // each block should have the address where the head of the list of "privatizable" blocks is kept
    // the only exception is a block for boot strap which is initialized when TLS is yet NULL
    nextPrivatizable = tlsBin? (Block*)(tlsBin + index) : NULL;
    TRACEF(( "[ScalableMalloc trace] Empty block %p is initialized, owner is %d, objectSize is %d, bumpPtr is %p\n",
             this, owner, objectSize, bumpPtr ));
}

Block *OrphanedBlocks::get(Bin* bin, unsigned int size)
{
    Block *result;
    MALLOC_ASSERT( bin, ASSERT_TEXT );
    unsigned int index = getIndex(size);
    result = bins[index].pop();
    if (result) {
        MALLOC_ITT_SYNC_ACQUIRED(bins+index);
        result->privatizeOrphaned(bin);
        STAT_increment(result->owner, index, allocBlockPublic);
    }
    return result;
}

void OrphanedBlocks::put(Bin* bin, Block *block)
{
    unsigned int index = getIndex(block->getSize());
    block->shareOrphaned(bin);
    MALLOC_ITT_SYNC_RELEASING(bins+index);
    bins[index].push(block);
}

void OrphanedBlocks::reset()
{
    for (uint32_t i=0; i<numBlockBinLimit; i++)
        new (bins+i) LifoList();
}

void FreeBlockPool::insertBlock(Block *block)
{
    size++;
    block->next = head;
    head = block;
    if (!tail)
        tail = block;
}

FreeBlockPool::ResOfGet FreeBlockPool::getBlock()
{
    Block *b = head;

    if (head) {
        size--;
        head = head->next;
        if (!head)
            tail = NULL;
        lastAccessMiss = false;
    } else
        lastAccessMiss = true;

    return ResOfGet(b, lastAccessMiss);
}

void FreeBlockPool::returnBlock(Block *block)
{
    MALLOC_ASSERT( size <= POOL_HIGH_MARK, ASSERT_TEXT );
    if (size == POOL_HIGH_MARK) {
        // release cold blocks and add hot one
        Block *headToFree = head,
              *helper;
        for (int i=0; i<POOL_LOW_MARK-2; i++)
            headToFree = headToFree->next;
        tail = headToFree;
        headToFree = headToFree->next;
        tail->next = NULL;
        size = POOL_LOW_MARK-1;
        // slab blocks from user pools not have valid backreference
        for (Block *currBl = headToFree; currBl; currBl = helper) {
            helper = currBl->next;
            if (!backend->inUserPool())
                removeBackRef(currBl->backRefIdx);
            backend->putSlabBlock(currBl);
        }
    }
    insertBlock(block);
}

bool FreeBlockPool::releaseAllBlocks()
{
    Block *helper;
    bool nonEmpty = size;

    for (Block *currBl = head; currBl; currBl=helper) {
        helper = currBl->next;
        // slab blocks in user's pools not have valid backRefIdx
        if (!backend->inUserPool())
            removeBackRef(currBl->backRefIdx);
        backend->putSlabBlock(currBl);
    }
    head = tail = NULL;
    size = 0;

    return nonEmpty;
}

/* We have a block give it back to the malloc block manager */
void Block::makeEmpty()
{
    // it is caller's responsibility to ensure no data is lost before calling this
    MALLOC_ASSERT( allocatedCount==0, ASSERT_TEXT );
    MALLOC_ASSERT( publicFreeList==NULL, ASSERT_TEXT );
    STAT_increment(owner, getIndex(objectSize), freeBlockBack);

    cleanBlockHeader();

    nextPrivatizable = NULL;

    objectSize = 0;
    owner.invalid();
    // for an empty block, bump pointer should point right after the end of the block
    bumpPtr = (FreeObject *)((uintptr_t)this + slabSize);
}

inline void Bin::setActiveBlock (Block *block)
{
//    MALLOC_ASSERT( bin, ASSERT_TEXT );
    MALLOC_ASSERT( block->owner.own(), ASSERT_TEXT );
    // it is the caller responsibility to keep bin consistence (i.e. ensure this block is in the bin list)
    activeBlk = block;
}

inline Block* Bin::setPreviousBlockActive()
{
    MALLOC_ASSERT( activeBlk, ASSERT_TEXT );
    Block* temp = activeBlk->previous;
    if( temp ) {
        MALLOC_ASSERT( temp->isFull == 0, ASSERT_TEXT );
        activeBlk = temp;
    }
    return temp;
}

FreeObject *Block::findObjectToFree(void *object) const
{
    FreeObject *objectToFree;
    // Due to aligned allocations, a pointer passed to scalable_free
    // might differ from the address of internally allocated object.
    // Small objects however should always be fine.
    if (objectSize <= maxSegregatedObjectSize)
        objectToFree = (FreeObject*)object;
    // "Fitting size" allocations are suspicious if aligned higher than naturally
    else {
        if ( ! isAligned(object,2*fittingAlignment) )
            // TODO: the above check is questionable - it gives false negatives in ~50% cases,
            //       so might even be slower in average than unconditional use of findAllocatedObject.
            // here it should be a "real" object
            objectToFree = (FreeObject*)object;
        else
            // here object can be an aligned address, so applying additional checks
            objectToFree = findAllocatedObject(object);
        MALLOC_ASSERT( isAligned(objectToFree,fittingAlignment), ASSERT_TEXT );
    }
    MALLOC_ASSERT( isProperlyPlaced(objectToFree), ASSERT_TEXT );

    return objectToFree;
}

void TLSData::release(MemoryPool *mPool)
{
    lloc.clean(&mPool->extMemPool);
    freeSlabBlocks.releaseAllBlocks();

    for (unsigned index = 0; index < numBlockBins; index++) {
        Block *activeBlk = bin[index].getActiveBlock();
        if (!activeBlk)
            continue;
        Block *threadlessBlock = activeBlk->previous;
        while (threadlessBlock) {
            Block *threadBlock = threadlessBlock->previous;
            if (threadlessBlock->empty()) {
                /* we destroy the thread, so not use its block pool */
                mPool->returnEmptyBlock(threadlessBlock, /*poolTheBlock=*/false);
            } else {
                mPool->orphanedBlocks.put(bin+index, threadlessBlock);
            }
            threadlessBlock = threadBlock;
        }
        threadlessBlock = activeBlk;
        while (threadlessBlock) {
            Block *threadBlock = threadlessBlock->next;
            if (threadlessBlock->empty()) {
                /* we destroy the thread, so not use its block pool */
                mPool->returnEmptyBlock(threadlessBlock, /*poolTheBlock=*/false);
            } else {
                mPool->orphanedBlocks.put(bin+index, threadlessBlock);
            }
            threadlessBlock = threadBlock;
        }
        bin[index].resetActiveBlock();
    }
}


#if MALLOC_CHECK_RECURSION
// TODO: Use deducated heap for this

/*
 * It's a special kind of allocation that can be used when malloc is
 * not available (either during startup or when malloc was already called and
 * we are, say, inside pthread_setspecific's call).
 * Block can contain objects of different sizes,
 * allocations are performed by moving bump pointer and increasing of object counter,
 * releasing is done via counter of objects allocated in the block
 * or moving bump pointer if releasing object is on a bound.
 */

class StartupBlock : public Block {
    size_t availableSize() const {
        return slabSize - ((uintptr_t)bumpPtr - (uintptr_t)this);
    }
    static StartupBlock *getBlock();
public:
    static FreeObject *allocate(size_t size);
    static size_t msize(void *ptr) { return *((size_t*)ptr - 1); }
    void free(void *ptr);
};

static MallocMutex startupMallocLock;
static StartupBlock *firstStartupBlock;

StartupBlock *StartupBlock::getBlock()
{
    BackRefIdx backRefIdx = BackRefIdx::newBackRef(/*largeObj=*/false);
    if (backRefIdx.isInvalid()) return NULL;

    StartupBlock *block = static_cast<StartupBlock*>(
        defaultMemPool->extMemPool.backend.getSlabBlock(1));
    if (!block) return NULL;

    block->cleanBlockHeader();
    setBackRef(backRefIdx, block);
    block->backRefIdx = backRefIdx;
    // use startupAllocObjSizeMark to mark objects from startup block marker
    block->objectSize = startupAllocObjSizeMark;
    block->bumpPtr = (FreeObject *)((uintptr_t)block + sizeof(StartupBlock));
    return block;
}

FreeObject *StartupBlock::allocate(size_t size)
{
    FreeObject *result;
    StartupBlock *newBlock = NULL;
    bool newBlockUnused = false;

    /* Objects must be aligned on their natural bounds,
       and objects bigger than word on word's bound. */
    size = alignUp(size, sizeof(size_t));
    // We need size of an object to implement msize.
    size_t reqSize = size + sizeof(size_t);
    // speculatively allocates newBlock to try avoid allocation while holding lock
    /* TODO: The function is called when malloc nested call is detected,
             so simultaneous usage from different threads seems unlikely.
             If pre-allocation is found useless, the code might be simplified. */
    if (!firstStartupBlock || firstStartupBlock->availableSize() < reqSize) {
        newBlock = StartupBlock::getBlock();
        if (!newBlock) return NULL;
    }
    {
        MallocMutex::scoped_lock scoped_cs(startupMallocLock);
        // Re-check whether we need a new block (conditions might have changed)
        if (!firstStartupBlock || firstStartupBlock->availableSize() < reqSize) {
            if (!newBlock) {
                newBlock = StartupBlock::getBlock();
                if (!newBlock) return NULL;
            }
            newBlock->next = (Block*)firstStartupBlock;
            if (firstStartupBlock)
                firstStartupBlock->previous = (Block*)newBlock;
            firstStartupBlock = newBlock;
        } else
            newBlockUnused = true;
        result = firstStartupBlock->bumpPtr;
        firstStartupBlock->allocatedCount++;
        firstStartupBlock->bumpPtr =
            (FreeObject *)((uintptr_t)firstStartupBlock->bumpPtr + reqSize);
    }
    if (newBlock && newBlockUnused)
        defaultMemPool->returnEmptyBlock(newBlock, /*poolTheBlock=*/false);

    // keep object size at the negative offset
    *((size_t*)result) = size;
    return (FreeObject*)((size_t*)result+1);
}

void StartupBlock::free(void *ptr)
{
    Block* blockToRelease = NULL;
    {
        MallocMutex::scoped_lock scoped_cs(startupMallocLock);

        MALLOC_ASSERT(firstStartupBlock, ASSERT_TEXT);
        MALLOC_ASSERT(startupAllocObjSizeMark==objectSize
                      && allocatedCount>0, ASSERT_TEXT);
        MALLOC_ASSERT((uintptr_t)ptr>=(uintptr_t)this+sizeof(StartupBlock)
                      && (uintptr_t)ptr+StartupBlock::msize(ptr)<=(uintptr_t)this+slabSize,
                      ASSERT_TEXT);
        if (0 == --allocatedCount) {
            if (this == firstStartupBlock)
                firstStartupBlock = (StartupBlock*)firstStartupBlock->next;
            if (previous)
                previous->next = next;
            if (next)
                next->previous = previous;
            blockToRelease = this;
        } else if ((uintptr_t)ptr + StartupBlock::msize(ptr) == (uintptr_t)bumpPtr) {
            // last object in the block released
            FreeObject *newBump = (FreeObject*)((size_t*)ptr - 1);
            MALLOC_ASSERT((uintptr_t)newBump>(uintptr_t)this+sizeof(StartupBlock),
                          ASSERT_TEXT);
            bumpPtr = newBump;
        }
    }
    if (blockToRelease) {
        blockToRelease->previous = blockToRelease->next = NULL;
        defaultMemPool->returnEmptyBlock(blockToRelease, /*poolTheBlock=*/false);
    }
}

#endif /* MALLOC_CHECK_RECURSION */

/********* End thread related code  *************/

/********* Library initialization *************/

//! Value indicating the state of initialization.
/* 0 = initialization not started.
 * 1 = initialization started but not finished.
 * 2 = initialization finished.
 * In theory, we only need values 0 and 2. But value 1 is nonetheless
 * useful for detecting errors in the double-check pattern.
 */
static intptr_t mallocInitialized;   // implicitly initialized to 0
static MallocMutex initMutex;

#include "../tbb/tbb_version.h"

/** The leading "\0" is here so that applying "strings" to the binary
    delivers a clean result. */
static char VersionString[] = "\0" TBBMALLOC_VERSION_STRINGS;

#if _XBOX || __TBB_WIN8UI_SUPPORT
bool GetBoolEnvironmentVariable(const char *) { return false; }
#else
bool GetBoolEnvironmentVariable(const char *name)
{
    if( const char* s = getenv(name) )
        return strcmp(s,"0") != 0;
    return false;
}
#endif

void AllocControlledMode::initReadEnv(const char *envName, intptr_t defaultVal)
{
    if (!setDone) {
#if !_XBOX && !__TBB_WIN8UI_SUPPORT
        const char *envVal = getenv(envName);
        if (envVal && !strcmp(envVal, "1"))
            val = 1;
        else
#endif
            val = defaultVal;
        setDone = true;
    }
}

void MemoryPool::initDefaultPool()
{
    long long hugePageSize = 0;
#if __linux__
    if (FILE *f = fopen("/proc/meminfo", "r")) {
        const int READ_BUF_SIZE = 100;
        char buf[READ_BUF_SIZE];
        MALLOC_ASSERT(sizeof(hugePageSize) >= 8,
                      "At least 64 bits required for keeping page size/numbers.");

        while (fgets(buf, READ_BUF_SIZE, f)) {
            if (1 == sscanf(buf, "Hugepagesize: %llu kB", &hugePageSize)) {
                hugePageSize *= 1024;
                break;
            }
        }
        fclose(f);
    }
#endif
    hugePages.init(hugePageSize);
}

inline bool isMallocInitialized() {
    // Load must have acquire fence; otherwise thread taking "initialized" path
    // might perform textually later loads *before* mallocInitialized becomes 2.
    return 2 == FencedLoad(mallocInitialized);
}

bool isMallocInitializedExt() {
    return isMallocInitialized();
}

/*
 * Allocator initialization routine;
 * it is called lazily on the very first scalable_malloc call.
 */
static void initMemoryManager()
{
    TRACEF(( "[ScalableMalloc trace] sizeof(Block) is %d (expected 128); sizeof(uintptr_t) is %d\n",
             sizeof(Block), sizeof(uintptr_t) ));
    MALLOC_ASSERT( 2*blockHeaderAlignment == sizeof(Block), ASSERT_TEXT );
    MALLOC_ASSERT( sizeof(FreeObject) == sizeof(void*), ASSERT_TEXT );

    bool initOk = defaultMemPool->
        extMemPool.init(0, NULL, NULL, scalableMallocPoolGranularity,
                        /*keepAllMemory=*/false, /*fixedPool=*/false);
// TODO: add error handling, and on error do something better than exit(1)
    if (!initOk || !initBackRefMaster(&defaultMemPool->extMemPool.backend)) {
        fprintf (stderr, "The memory manager cannot access sufficient memory to initialize; exiting \n");
        exit(1);
    }
    ThreadId::init();      // Create keys for thread id
    MemoryPool::initDefaultPool();
#if COLLECT_STATISTICS
    initStatisticsCollection();
#endif
}

//! Ensures that initMemoryManager() is called once and only once.
/** Does not return until initMemoryManager() has been completed by a thread.
    There is no need to call this routine if mallocInitialized==2 . */
static void doInitialization()
{
    MallocMutex::scoped_lock lock( initMutex );
    if (mallocInitialized!=2) {
        MALLOC_ASSERT( mallocInitialized==0, ASSERT_TEXT );
        mallocInitialized = 1;
        RecursiveMallocCallProtector scoped;
        initMemoryManager();
#ifdef  MALLOC_EXTRA_INITIALIZATION
        MALLOC_EXTRA_INITIALIZATION;
#endif
#if MALLOC_CHECK_RECURSION
        RecursiveMallocCallProtector::detectNaiveOverload();
#endif
        MALLOC_ASSERT( mallocInitialized==1, ASSERT_TEXT );
        // Store must have release fence, otherwise mallocInitialized==2
        // might become remotely visible before side effects of
        // initMemoryManager() become remotely visible.
        FencedStore( mallocInitialized, 2 );
        if( GetBoolEnvironmentVariable("TBB_VERSION") ) {
            fputs(VersionString+1,stderr);
            hugePages.printStatus();
        }
    }
    /* It can't be 0 or I would have initialized it */
    MALLOC_ASSERT( mallocInitialized==2, ASSERT_TEXT );
}

/********* End library initialization *************/

/********* The malloc show begins     *************/


FreeObject *Block::allocateFromFreeList()
{
    FreeObject *result;

    if (!freeList) return NULL;

    result = freeList;
    MALLOC_ASSERT( result, ASSERT_TEXT );

    freeList = result->next;
    MALLOC_ASSERT( allocatedCount < (slabSize-sizeof(Block))/objectSize, ASSERT_TEXT );
    allocatedCount++;
    STAT_increment(owner, getIndex(objectSize), allocFreeListUsed);

    return result;
}

FreeObject *Block::allocateFromBumpPtr()
{
    FreeObject *result = bumpPtr;
    if (result) {
        bumpPtr = (FreeObject *) ((uintptr_t) bumpPtr - objectSize);
        if ( (uintptr_t)bumpPtr < (uintptr_t)this+sizeof(Block) ) {
            bumpPtr = NULL;
        }
        MALLOC_ASSERT( allocatedCount < (slabSize-sizeof(Block))/objectSize, ASSERT_TEXT );
        allocatedCount++;
        STAT_increment(owner, getIndex(objectSize), allocBumpPtrUsed);
    }
    return result;
}

inline FreeObject* Block::allocate()
{
    MALLOC_ASSERT( owner.own(), ASSERT_TEXT );

    /* for better cache locality, first looking in the free list. */
    if ( FreeObject *result = allocateFromFreeList() ) {
        return result;
    }
    MALLOC_ASSERT( !freeList, ASSERT_TEXT );

    /* if free list is empty, try thread local bump pointer allocation. */
    if ( FreeObject *result = allocateFromBumpPtr() ) {
        return result;
    }
    MALLOC_ASSERT( !bumpPtr, ASSERT_TEXT );

    /* the block is considered full. */
    isFull = 1;
    return NULL;
}

void Bin::moveBlockToBinFront(Block *block)
{
    /* move the block to the front of the bin */
    if (block == activeBlk) return;
    outofTLSBin(block);
    pushTLSBin(block);
}

void Bin::processLessUsedBlock(MemoryPool *memPool, Block *block)
{
    if (block != activeBlk) {
        /* We are not actively using this block; return it to the general block pool */
        outofTLSBin(block);
        memPool->returnEmptyBlock(block, /*poolTheBlock=*/true);
    } else {
        /* all objects are free - let's restore the bump pointer */
        block->restoreBumpPtr();
    }
}

template<int LOW_MARK, int HIGH_MARK>
bool LocalLOC<LOW_MARK, HIGH_MARK>::put(LargeMemoryBlock *object, ExtMemoryPool *extMemPool)
{
    const size_t size = object->unalignedSize;
    if (size > MAX_TOTAL_SIZE)
        return false;

    totalSize += size;
    object->prev = NULL;
    object->next = head;
    if (head) head->prev = object;
    head = object;
    if (!tail) tail = object;
    numOfBlocks++;
    MALLOC_ASSERT(!tail->next, ASSERT_TEXT);
    // must meet both size and number of cached objects constrains
    if (totalSize > MAX_TOTAL_SIZE || numOfBlocks >= HIGH_MARK) {
        // scanning from tail until meet conditions
        while (totalSize > MAX_TOTAL_SIZE || numOfBlocks > LOW_MARK) {
            totalSize -= tail->unalignedSize;
            numOfBlocks--;
            tail = tail->prev;
        }
        LargeMemoryBlock *headToRelease = tail->next;
        tail->next = NULL;

        extMemPool->freeLargeObjectList(headToRelease);
    }
    lastUsedOSCallsCnt = lastSeenOSCallsCnt;
    return true;
}

template<int LOW_MARK, int HIGH_MARK>
LargeMemoryBlock *LocalLOC<LOW_MARK, HIGH_MARK>::get(size_t size)
{
    if (lastUsedOSCallsCnt != lastSeenOSCallsCnt)
        lastUsedOSCallsCnt = lastSeenOSCallsCnt;

    for (LargeMemoryBlock *curr = head; curr; curr=curr->next) {
        if (curr->unalignedSize == size) {
            LargeMemoryBlock *res = curr;
            if (curr->next)
                curr->next->prev = curr->prev;
            else
                tail = curr->prev;
            if (curr->prev)
                curr->prev->next = curr->next;
            else
                head = curr->next;
            totalSize -= size;
            numOfBlocks--;
            return res;
        }
    }
    return NULL;
}

template<int LOW_MARK, int HIGH_MARK>
bool LocalLOC<LOW_MARK, HIGH_MARK>::clean(ExtMemoryPool *extMemPool)
{
    bool released = numOfBlocks;

    if (numOfBlocks)
        extMemPool->freeLargeObjectList(head);
    head = tail = NULL;
    numOfBlocks = 0;
    totalSize = 0;
    return released;
}

template<int LOW_MARK, int HIGH_MARK>
void LocalLOC<LOW_MARK, HIGH_MARK>::allocatorCalledHook(ExtMemoryPool *extMemPool)
{
    intptr_t currCnt = extMemPool->backend.askMemFromOSCounter.get();

    // clean the cache iff there was OS memory request since last hook call
    // and the cache was not touched since previous OS memory request
    if (currCnt != lastSeenOSCallsCnt && lastUsedOSCallsCnt != lastSeenOSCallsCnt
        && head)
        clean(extMemPool);
    lastSeenOSCallsCnt = currCnt;
}

void *MemoryPool::getFromLLOCache(TLSData* tls, size_t size, size_t alignment)
{
    LargeMemoryBlock *lmb = NULL;

    size_t headersSize = sizeof(LargeMemoryBlock)+sizeof(LargeObjectHdr);
    size_t allocationSize = LargeObjectCache::alignToBin(size+headersSize+alignment);
    if (allocationSize < size) // allocationSize is wrapped around after alignToBin
        return NULL;

    if (tls)
        lmb = tls->lloc.get(allocationSize);
    if (!lmb)
        lmb = extMemPool.mallocLargeObject(allocationSize);

    if (lmb) {
        void *alignedArea = (void*)alignUp((uintptr_t)lmb+headersSize, alignment);
        LargeObjectHdr *header = (LargeObjectHdr*)alignedArea-1;
        header->memoryBlock = lmb;
        header->backRefIdx = lmb->backRefIdx;
        setBackRef(header->backRefIdx, header);

        lmb->objectSize = size;

        MALLOC_ASSERT( isLargeObject(alignedArea), ASSERT_TEXT );

        return alignedArea;
    }
    return NULL;
}

void MemoryPool::putToLLOCache(TLSData *tls, void *object)
{
    LargeObjectHdr *header = (LargeObjectHdr*)object - 1;
    // overwrite backRefIdx to simplify double free detection
    header->backRefIdx = BackRefIdx();

    if (!tls || !tls->lloc.put(header->memoryBlock, &extMemPool))
        extMemPool.freeLargeObject(header->memoryBlock);
}

// called on each allocator call
void MemoryPool::allocatorCalledHook(TLSData *tls)
{
    // TODO: clean freeSlabBlocks as well
    tls->lloc.allocatorCalledHook(&extMemPool);
}

#if USE_PTHREAD && (__TBB_SOURCE_DIRECTLY_INCLUDED || __TBB_USE_DLOPEN_REENTRANCY_WORKAROUND)

/* Decrease race interval between dynamic library unloading and pthread key
   destructor. Protect only Pthreads with supported unloading. */
class ShutdownSync {
/* flag is the number of threads in pthread key dtor body
   (i.e., between threadDtorStart() and threadDtorDone())
   or the signal to skip dtor, if flag < 0 */
    intptr_t flag;
    static const intptr_t skipDtor = INTPTR_MIN/2;
public:
/* Suppose that 2*abs(skipDtor) or more threads never call threadExitStart()
   simultaneously, so flag is never becomes negative because of that. */
    bool threadDtorStart() {
        if (flag < 0)
            return false;
        if (AtomicIncrement(flag) <= 0) { // note that new value returned
            AtomicAdd(flag, -1);  // flag is spoiled by us, restore it
            return false;
        }
        return true;
    }
    void threadDtorDone() {
        AtomicAdd(flag, -1);
    }
    void processExit() {
        if (AtomicAdd(flag, skipDtor) != 0)
            SpinWaitUntilEq(flag, skipDtor);
    }
};

#else

class ShutdownSync {
public:
    bool threadDtorStart() { return true; }
    void threadDtorDone() { }
    void processExit() { }
};

#endif // USE_PTHREAD && (__TBB_SOURCE_DIRECTLY_INCLUDED || __TBB_USE_DLOPEN_REENTRANCY_WORKAROUND)

static ShutdownSync shutdownSync;

/*
 * All aligned allocations fall into one of the following categories:
 *  1. if both request size and alignment are <= maxSegregatedObjectSize,
 *       we just align the size up, and request this amount, because for every size
 *       aligned to some power of 2, the allocated object is at least that aligned.
 * 2. for size<minLargeObjectSize, check if already guaranteed fittingAlignment is enough.
 * 3. if size+alignment<minLargeObjectSize, we take an object of fittingSizeN and align
 *       its address up; given such pointer, scalable_free could find the real object.
 *       Wrapping of size+alignment is impossible because maximal allowed
 *       alignment plus minLargeObjectSize can't lead to wrapping.
 * 4. otherwise, aligned large object is allocated.
 */
static void *allocateAligned(MemoryPool *memPool, size_t size, size_t alignment)
{
    MALLOC_ASSERT( isPowerOfTwo(alignment), ASSERT_TEXT );

    if (!isMallocInitialized()) doInitialization();

    void *result;
    if (size<=maxSegregatedObjectSize && alignment<=maxSegregatedObjectSize)
        result = internalPoolMalloc(memPool, alignUp(size? size: sizeof(size_t), alignment));
    else if (size<minLargeObjectSize) {
        if (alignment<=fittingAlignment)
            result = internalPoolMalloc(memPool, size);
        else if (size+alignment < minLargeObjectSize) {
            void *unaligned = internalPoolMalloc(memPool, size+alignment);
            if (!unaligned) return NULL;
            result = alignUp(unaligned, alignment);
        } else
            goto LargeObjAlloc;
    } else {
    LargeObjAlloc:
        /* This can be the first allocation call. */
        if (!isMallocInitialized())
            doInitialization();
        TLSData *tls = memPool->getTLS(/*create=*/true);
        memPool->allocatorCalledHook(tls);
        // take into account only alignment that are higher then natural
        result =
            memPool->getFromLLOCache(tls, size, largeObjectAlignment>alignment?
                                               largeObjectAlignment: alignment);
    }

    MALLOC_ASSERT( isAligned(result, alignment), ASSERT_TEXT );
    return result;
}

static void *reallocAligned(MemoryPool *memPool, void *ptr,
                            size_t size, size_t alignment = 0)
{
    void *result;
    size_t copySize;

    if (isLargeObject(ptr)) {
        LargeMemoryBlock* lmb = ((LargeObjectHdr *)ptr - 1)->memoryBlock;
        copySize = lmb->unalignedSize-((uintptr_t)ptr-(uintptr_t)lmb);
        if (size <= copySize && (0==alignment || isAligned(ptr, alignment))) {
            lmb->objectSize = size;
            return ptr;
        } else {
            copySize = lmb->objectSize;
            result = alignment ? allocateAligned(memPool, size, alignment) :
                internalPoolMalloc(memPool, size);
        }
    } else {
        Block* block = (Block *)alignDown(ptr, slabSize);
        copySize = block->getSize();
        if (size <= copySize && (0==alignment || isAligned(ptr, alignment))) {
            return ptr;
        } else {
            result = alignment ? allocateAligned(memPool, size, alignment) :
                internalPoolMalloc(memPool, size);
        }
    }
    if (result) {
        memcpy(result, ptr, copySize<size? copySize: size);
        internalPoolFree(memPool, ptr);
    }
    return result;
}

/* A predicate checks if an object is properly placed inside its block */
inline bool Block::isProperlyPlaced(const void *object) const
{
    return 0 == ((uintptr_t)this + slabSize - (uintptr_t)object) % objectSize;
}

/* Finds the real object inside the block */
FreeObject *Block::findAllocatedObject(const void *address) const
{
    // calculate offset from the end of the block space
    uint16_t offset = (uintptr_t)this + slabSize - (uintptr_t)address;
    MALLOC_ASSERT( offset<=slabSize-sizeof(Block), ASSERT_TEXT );
    // find offset difference from a multiple of allocation size
    offset %= objectSize;
    // and move the address down to where the real object starts.
    return (FreeObject*)((uintptr_t)address - (offset? objectSize-offset: 0));
}

/*
 * Bad dereference caused by a foreign pointer is possible only here, not earlier in call chain.
 * Separate function isolates SEH code, as it has bad influence on compiler optimization.
 */
static inline BackRefIdx safer_dereference (const BackRefIdx *ptr)
{
    BackRefIdx id;
#if _MSC_VER
    __try {
#endif
        id = *ptr;
#if _MSC_VER
    } __except( GetExceptionCode() == EXCEPTION_ACCESS_VIOLATION?
                EXCEPTION_EXECUTE_HANDLER : EXCEPTION_CONTINUE_SEARCH ) {
        id = BackRefIdx();
    }
#endif
    return id;
}

bool isLargeObject(void *object)
{
    if (!isAligned(object, largeObjectAlignment))
        return false;
    LargeObjectHdr *header = (LargeObjectHdr*)object - 1;
    BackRefIdx idx = safer_dereference(&header->backRefIdx);

    return idx.isLargeObject()
        // in valid LargeObjectHdr memoryBlock points somewhere before header
        // TODO: more strict check
        && (uintptr_t)header->memoryBlock < (uintptr_t)header
        && getBackRef(idx) == header;
}

static inline bool isSmallObject (void *ptr)
{
    void* expected = alignDown(ptr, slabSize);
    const BackRefIdx* idx = ((Block*)expected)->getBackRef();

    return expected == getBackRef(safer_dereference(idx));
}

/**** Check if an object was allocated by scalable_malloc ****/
static inline bool isRecognized (void* ptr)
{
    return isLargeObject(ptr) || isSmallObject(ptr);
}

static inline void freeSmallObject(MemoryPool *memPool, TLSData *tls, void *object)
{
    /* mask low bits to get the block */
    Block *block = (Block *)alignDown(object, slabSize);
    MALLOC_ASSERT( block->checkFreePrecond(object),
                   "Possible double free or heap corruption." );

#if MALLOC_CHECK_RECURSION
    if (block->isStartupAllocObject()) {
        ((StartupBlock *)block)->free(object);
        return;
    }
#endif
    if (block->ownBlock())
        block->freeOwnObject(memPool, tls, object);
    else { /* Slower path to add to the shared list, the allocatedCount is updated by the owner thread in malloc. */
        FreeObject *objectToFree = block->findObjectToFree(object);
        block->freePublicObject(objectToFree);
    }
}

static void *internalPoolMalloc(MemoryPool* memPool, size_t size)
{
    Bin* bin;
    Block * mallocBlock;

    if (!memPool) return NULL;

    if (!size) size = sizeof(size_t);

    TLSData *tls = memPool->getTLS(/*create=*/true);
    memPool->allocatorCalledHook(tls);
    /*
     * Use Large Object Allocation
     */
    if (size >= minLargeObjectSize)
        return memPool->getFromLLOCache(tls, size, largeObjectAlignment);

    /*
     * Get an element in thread-local array corresponding to the given size;
     * It keeps ptr to the active block for allocations of this size
     */
    bin = memPool->getAllocationBin(tls, size);
    if ( !bin ) return NULL;

    /* Get a block to try to allocate in. */
    for( mallocBlock = bin->getActiveBlock(); mallocBlock;
         mallocBlock = bin->setPreviousBlockActive() ) // the previous block should be empty enough
    {
        if( FreeObject *result = mallocBlock->allocate() )
            return result;
    }

    /*
     * else privatize publicly freed objects in some block and allocate from it
     */
    mallocBlock = bin->getPublicFreeListBlock();
    if (mallocBlock) {
        if (mallocBlock->emptyEnoughToUse()) {
            bin->moveBlockToBinFront(mallocBlock);
        }
        MALLOC_ASSERT( mallocBlock->freeListNonNull(), ASSERT_TEXT );
        if ( FreeObject *result = mallocBlock->allocateFromFreeList() )
            return result;
        /* Else something strange happened, need to retry from the beginning; */
        TRACEF(( "[ScalableMalloc trace] Something is wrong: no objects in public free list; reentering.\n" ));
        return internalPoolMalloc(memPool, size);
    }

    /*
     * no suitable own blocks, try to get a partial block that some other thread has discarded.
     */
    mallocBlock = memPool->orphanedBlocks.get(bin, size);
    while (mallocBlock) {
        bin->pushTLSBin(mallocBlock);
        bin->setActiveBlock(mallocBlock); // TODO: move under the below condition?
        if( FreeObject *result = mallocBlock->allocate() )
            return result;
        mallocBlock = memPool->orphanedBlocks.get(bin, size);
    }

    /*
     * else try to get a new empty block
     */
    mallocBlock = memPool->getEmptyBlock(size);
    if (mallocBlock) {
        bin->pushTLSBin(mallocBlock);
        bin->setActiveBlock(mallocBlock);
        if( FreeObject *result = mallocBlock->allocate() )
            return result;
        /* Else something strange happened, need to retry from the beginning; */
        TRACEF(( "[ScalableMalloc trace] Something is wrong: no objects in empty block; reentering.\n" ));
        return internalPoolMalloc(memPool, size);
    }
    /*
     * else nothing works so return NULL
     */
    TRACEF(( "[ScalableMalloc trace] No memory found, returning NULL.\n" ));
    return NULL;
}

static bool internalPoolFree(MemoryPool *memPool, void *object)
{
    if (!memPool || !object) return false;

    // The library is initialized at allocation call, so releasing while
    // not initialized means foreign object is releasing.
    MALLOC_ASSERT(isMallocInitialized(), ASSERT_TEXT);
    MALLOC_ASSERT(memPool->extMemPool.userPool() || isRecognized(object),
                  "Invalid pointer in pool_free detected.");
    TLSData *tls = memPool->getTLS(/*create=*/false);
    if (tls) memPool->allocatorCalledHook(tls);

    if (isLargeObject(object))
        memPool->putToLLOCache(tls, object);
    else
        freeSmallObject(memPool, tls, object);
    return true;
}

static void *internalMalloc(size_t size)
{
    if (!size) size = sizeof(size_t);

#if MALLOC_CHECK_RECURSION
    if (RecursiveMallocCallProtector::sameThreadActive())
        return size<minLargeObjectSize? StartupBlock::allocate(size) :
            // nested allocation, so skip tls
            (FreeObject*)defaultMemPool->getFromLLOCache(NULL, size, slabSize);
#endif

    if (!isMallocInitialized())
        doInitialization();

    return internalPoolMalloc(defaultMemPool, size);
}

static void internalFree(void *object)
{
    internalPoolFree(defaultMemPool, object);
}

static size_t internalMsize(void* ptr)
{
    if (ptr) {
        MALLOC_ASSERT(isRecognized(ptr), "Invalid pointer in scalable_msize detected.");
        if (isLargeObject(ptr)) {
            LargeMemoryBlock* lmb = ((LargeObjectHdr*)ptr - 1)->memoryBlock;
            return lmb->objectSize;
        } else {
            Block* block = (Block *)alignDown(ptr, slabSize);
#if MALLOC_CHECK_RECURSION
            size_t size = block->getSize()? block->getSize() : StartupBlock::msize(ptr);
#else
            size_t size = block->getSize();
#endif
            MALLOC_ASSERT(size>0 && size<minLargeObjectSize, ASSERT_TEXT);
            return size;
        }
    }
    errno = EINVAL;
    // Unlike _msize, return 0 in case of parameter error.
    // Returning size_t(-1) looks more like the way to troubles.
    return 0;
}

} // namespace internal

using namespace rml::internal;

// legacy entry point saved for compatibility with binaries complied
// with pre-6003 versions of TBB
rml::MemoryPool *pool_create(intptr_t pool_id, const MemPoolPolicy *policy)
{
    rml::MemoryPool *pool;
    MemPoolPolicy pol(policy->pAlloc, policy->pFree, policy->granularity);

    pool_create_v1(pool_id, &pol, &pool);
    return pool;
}

rml::MemPoolError pool_create_v1(intptr_t pool_id, const MemPoolPolicy *policy,
                                 rml::MemoryPool **pool)
{
    if ( !policy->pAlloc || policy->version<MemPoolPolicy::TBBMALLOC_POOL_VERSION
         // empty pFree allowed only for fixed pools
         || !(policy->fixedPool || policy->pFree) ) {
        *pool = NULL;
        return INVALID_POLICY;
    }
    if ( policy->version>MemPoolPolicy::TBBMALLOC_POOL_VERSION // future versions are not supported
         // new flags can be added in place of reserved, but default
         // behaviour must be supported by this version
         || policy->reserved ) {
        *pool = NULL;
        return UNSUPPORTED_POLICY;
    }
    if (!isMallocInitialized())
        doInitialization();

    rml::internal::MemoryPool *memPool =
        (rml::internal::MemoryPool*)internalMalloc((sizeof(rml::internal::MemoryPool)));
    if (!memPool) {
        *pool = NULL;
        return NO_MEMORY;
    }
    memset(memPool, 0, sizeof(rml::internal::MemoryPool));
    if (!memPool->init(pool_id, policy)) {
        internalFree(memPool);
        *pool = NULL;
        return NO_MEMORY;
    }

    *pool = (rml::MemoryPool*)memPool;
    return POOL_OK;
}

bool pool_destroy(rml::MemoryPool* memPool)
{
    if (!memPool) return false;
    ((rml::internal::MemoryPool*)memPool)->destroy();
    internalFree(memPool);

    return true;
}

bool pool_reset(rml::MemoryPool* memPool)
{
    if (!memPool) return false;

    ((rml::internal::MemoryPool*)memPool)->reset();
    return true;
}

void *pool_malloc(rml::MemoryPool* mPool, size_t size)
{
    return internalPoolMalloc((rml::internal::MemoryPool*)mPool, size);
}

void *pool_realloc(rml::MemoryPool* mPool, void *object, size_t size)
{
    if (!object)
        return internalPoolMalloc((rml::internal::MemoryPool*)mPool, size);
    if (!size) {
        internalPoolFree((rml::internal::MemoryPool*)mPool, object);
        return NULL;
    }
    return reallocAligned((rml::internal::MemoryPool*)mPool, object, size, 0);
}

void *pool_aligned_malloc(rml::MemoryPool* mPool, size_t size, size_t alignment)
{
    if (!isPowerOfTwo(alignment) || 0==size)
        return NULL;

    return allocateAligned((rml::internal::MemoryPool*)mPool, size, alignment);
}

void *pool_aligned_realloc(rml::MemoryPool* memPool, void *ptr, size_t size, size_t alignment)
{
    if (!isPowerOfTwo(alignment))
        return NULL;
    rml::internal::MemoryPool *mPool = (rml::internal::MemoryPool*)memPool;
    void *tmp;

    if (!ptr)
        tmp = allocateAligned(mPool, size, alignment);
    else if (!size) {
        internalPoolFree(mPool, ptr);
        return NULL;
    } else
        tmp = reallocAligned(mPool, ptr, size, alignment);

    return tmp;
}

bool pool_free(rml::MemoryPool *mPool, void *object)
{
    return internalPoolFree((rml::internal::MemoryPool*)mPool, object);
}

} // namespace rml

using namespace rml::internal;

#if MALLOC_TRACE
static unsigned int threadGoingDownCount = 0;
#endif

/*
 * When a thread is shutting down this routine should be called to remove all the thread ids
 * from the malloc blocks and replace them with a NULL thread id.
 *
 * For pthreads, the function is set as a callback in pthread_key_create for TLS bin.
 * For non-NULL keys it will be automatically called at thread exit with the key value
 * as the argument.
 *
 * for Windows, it should be called directly e.g. from DllMain
*/
void mallocThreadShutdownNotification(void* arg)
{
    // Check whether TLS has been initialized
    if (!isMallocInitialized()) return;

    TRACEF(( "[ScalableMalloc trace] Thread id %d blocks return start %d\n",
             getThreadId(),  threadGoingDownCount++ ));
#if USE_WINTHREAD
    suppress_unused_warning(arg);
    MallocMutex::scoped_lock lock(MemoryPool::memPoolListLock);
    // The routine is called once per thread, need to walk through all pools on Windows
    for (MemoryPool *memPool = defaultMemPool; memPool; memPool = memPool->next)
        if (TLSData *tls = memPool->getTLS(/*create=*/false))
            memPool->processThreadShutdown(tls);
#else
    if (!shutdownSync.threadDtorStart()) return;
    // The routine is called for each memPool, just need to get memPool from TLSData.
    TLSData *tls = (TLSData*)arg;
    tls->getMemPool()->processThreadShutdown(tls);
    shutdownSync.threadDtorDone();
#endif

    TRACEF(( "[ScalableMalloc trace] Thread id %d blocks return end\n", getThreadId() ));
}

#if USE_WINTHREAD
extern "C" void __TBB_mallocThreadShutdownNotification()
{
    mallocThreadShutdownNotification(NULL);
}
#endif

extern "C" void __TBB_mallocProcessShutdownNotification()
{
    if (!isMallocInitialized()) return;

#if __TBB_MALLOC_LOCACHE_STAT
    printf("cache hit ratio %f, size hit %f\n",
           1.*cacheHits/mallocCalls, 1.*memHitKB/memAllocKB);
    defaultMemPool->extMemPool.loc.reportStat(stdout);
#endif
    shutdownSync.processExit();
#if __TBB_SOURCE_DIRECTLY_INCLUDED
/* Pthread keys must be deleted as soon as possible to not call key dtor
   on thread termination when then the tbbmalloc code can be already unloaded.
*/
    defaultMemPool->destroy();
    destroyBackRefMaster(&defaultMemPool->extMemPool.backend);
    ThreadId::destroy();      // Delete key for thread id
#elif __TBB_USE_DLOPEN_REENTRANCY_WORKAROUND
/* In most cases we prevent unloading tbbmalloc, and don't clean up memory
   on process shutdown. When impossible to prevent, library unload results
   in shutdown notification, and it makes sense to release unused memory
   at that point (we can't release all memory because it's possible that
   it will be accessed after this point).
   TODO: better support systems where we can't prevent unloading by removing
   pthread destructors and releasing caches.
 */
    defaultMemPool->extMemPool.hardCachesCleanup();
#endif // __TBB_SOURCE_DIRECTLY_INCLUDED

#if COLLECT_STATISTICS
    ThreadId nThreads = ThreadIdCount;
    for( int i=1; i<=nThreads && i<MAX_THREADS; ++i )
        STAT_print(i);
#endif
}

extern "C" void * scalable_malloc(size_t size)
{
    void *ptr = internalMalloc(size);
    if (!ptr) errno = ENOMEM;
    return ptr;
}

extern "C" void scalable_free (void *object) {
    internalFree(object);
}

/*
 * A variant that provides additional memory safety, by checking whether the given address
 * was obtained with this allocator, and if not redirecting to the provided alternative call.
 */
extern "C" void safer_scalable_free (void *object, void (*original_free)(void*))
{
    if (!object)
        return;

    // must check 1st for large object, because small object check touches 4 pages on left,
    // and it can be unaccessable
    if (isLargeObject(object)) {
        TLSData *tls = defaultMemPool->getTLS(/*create=*/false);
        if (tls) defaultMemPool->allocatorCalledHook(tls);

        defaultMemPool->putToLLOCache(tls, object);
    } else if (isSmallObject(object)) {
        TLSData *tls = defaultMemPool->getTLS(/*create=*/false);
        if (tls) defaultMemPool->allocatorCalledHook(tls);

        freeSmallObject(defaultMemPool, tls, object);
    } else if (original_free)
        original_free(object);
}

/********* End the free code        *************/

/********* Code for scalable_realloc       ***********/

/*
 * From K&R
 * "realloc changes the size of the object pointed to by p to size. The contents will
 * be unchanged up to the minimum of the old and the new sizes. If the new size is larger,
 * the new space is uninitialized. realloc returns a pointer to the new space, or
 * NULL if the request cannot be satisfied, in which case *p is unchanged."
 *
 */
extern "C" void* scalable_realloc(void* ptr, size_t size)
{
    void *tmp;

    if (!ptr)
        tmp = internalMalloc(size);
    else if (!size) {
        internalFree(ptr);
        return NULL;
    } else
        tmp = reallocAligned(defaultMemPool, ptr, size, 0);

    if (!tmp) errno = ENOMEM;
    return tmp;
}

/*
 * A variant that provides additional memory safety, by checking whether the given address
 * was obtained with this allocator, and if not redirecting to the provided alternative call.
 */
extern "C" void* safer_scalable_realloc (void* ptr, size_t sz, void* original_realloc)
{
    void *tmp; // TODO: fix warnings about uninitialized use of tmp

    if (!ptr) {
        tmp = internalMalloc(sz);
    } else if (isRecognized(ptr)) {
        if (!sz) {
            internalFree(ptr);
            return NULL;
        } else {
            tmp = reallocAligned(defaultMemPool, ptr, sz, 0);
        }
    }
#if USE_WINTHREAD
    else if (original_realloc && sz) {
        orig_ptrs *original_ptrs = static_cast<orig_ptrs*>(original_realloc);
        if ( original_ptrs->orig_msize ){
            size_t oldSize = original_ptrs->orig_msize(ptr);
            tmp = internalMalloc(sz);
            if (tmp) {
                memcpy(tmp, ptr, sz<oldSize? sz : oldSize);
                if ( original_ptrs->orig_free ){
                    original_ptrs->orig_free( ptr );
                }
            }
        } else
            tmp = NULL;
    }
#else
    else if (original_realloc) {
        typedef void* (*realloc_ptr_t)(void*,size_t);
        realloc_ptr_t original_realloc_ptr;
        (void *&)original_realloc_ptr = original_realloc;
        tmp = original_realloc_ptr(ptr,sz);
    }
#endif
    else tmp = NULL;

    if (!tmp) errno = ENOMEM;
    return tmp;
}

/********* End code for scalable_realloc   ***********/

/********* Code for scalable_calloc   ***********/

/*
 * From K&R
 * calloc returns a pointer to space for an array of nobj objects,
 * each of size size, or NULL if the request cannot be satisfied.
 * The space is initialized to zero bytes.
 *
 */

extern "C" void * scalable_calloc(size_t nobj, size_t size)
{
    size_t arraySize = nobj * size;
    void* result = internalMalloc(arraySize);
    if (result)
        memset(result, 0, arraySize);
    else
        errno = ENOMEM;
    return result;
}

/********* End code for scalable_calloc   ***********/

/********* Code for aligned allocation API **********/

extern "C" int scalable_posix_memalign(void **memptr, size_t alignment, size_t size)
{
    if ( !isPowerOfTwoMultiple(alignment, sizeof(void*)) )
        return EINVAL;
    void *result = allocateAligned(defaultMemPool, size, alignment);
    if (!result)
        return ENOMEM;
    *memptr = result;
    return 0;
}

extern "C" void * scalable_aligned_malloc(size_t size, size_t alignment)
{
    if (!isPowerOfTwo(alignment) || 0==size) {
        errno = EINVAL;
        return NULL;
    }
    void *tmp = allocateAligned(defaultMemPool, size, alignment);
    if (!tmp) errno = ENOMEM;
    return tmp;
}

extern "C" void * scalable_aligned_realloc(void *ptr, size_t size, size_t alignment)
{
    if (!isPowerOfTwo(alignment)) {
        errno = EINVAL;
        return NULL;
    }
    void *tmp;

    if (!ptr)
        tmp = allocateAligned(defaultMemPool, size, alignment);
    else if (!size) {
        scalable_free(ptr);
        return NULL;
    } else
        tmp = reallocAligned(defaultMemPool, ptr, size, alignment);

    if (!tmp) errno = ENOMEM;
    return tmp;
}

extern "C" void * safer_scalable_aligned_realloc(void *ptr, size_t size, size_t alignment, void* orig_function)
{
    /* corner cases left out of reallocAligned to not deal with errno there */
    if (!isPowerOfTwo(alignment)) {
        errno = EINVAL;
        return NULL;
    }
    void *tmp = NULL;

    if (!ptr) {
        tmp = allocateAligned(defaultMemPool, size, alignment);
    } else if (isRecognized(ptr)) {
        if (!size) {
            internalFree(ptr);
            return NULL;
        } else {
            tmp = reallocAligned(defaultMemPool, ptr, size, alignment);
        }
    }
#if USE_WINTHREAD
    else {
        orig_ptrs *original_ptrs = static_cast<orig_ptrs*>(orig_function);
        if (size) {
            // Without orig_msize, we can't do anything with this.
            // Just keeping old pointer.
            if ( original_ptrs->orig_msize ){
                size_t oldSize = original_ptrs->orig_msize(ptr);
                tmp = allocateAligned(defaultMemPool, size, alignment);
                if (tmp) {
                    memcpy(tmp, ptr, size<oldSize? size : oldSize);
                    if ( original_ptrs->orig_free ){
                        original_ptrs->orig_free( ptr );
                    }
                }
            }
        } else {
            if ( original_ptrs->orig_free ){
                original_ptrs->orig_free( ptr );
            }
            return NULL;
        }
    }
#else
    // As original_realloc can't align result, and there is no way to find
    // size of reallocating object, we are giving up.
    suppress_unused_warning(orig_function);
#endif
    if (!tmp) errno = ENOMEM;
    return tmp;
}

extern "C" void scalable_aligned_free(void *ptr)
{
    internalFree(ptr);
}

/********* end code for aligned allocation API **********/

/********* Code for scalable_msize       ***********/

/*
 * Returns the size of a memory block allocated in the heap.
 */
extern "C" size_t scalable_msize(void* ptr)
{
    return internalMsize(ptr);
}

/*
 * A variant that provides additional memory safety, by checking whether the given address
 * was obtained with this allocator, and if not redirecting to the provided alternative call.
 */
extern "C" size_t safer_scalable_msize (void *object, size_t (*original_msize)(void*))
{
    if (object) {
        // Check if the memory was allocated by scalable_malloc
        if (isRecognized(object))
            return internalMsize(object);
        else if (original_msize)
            return original_msize(object);
    }
    // object is NULL or unknown
    errno = EINVAL;
    return 0;
}

/********* End code for scalable_msize   ***********/

extern "C" int scalable_allocation_mode(int param, intptr_t value)
{
#if __linux__
    if (param == USE_HUGE_PAGES)
        switch (value) {
        case 0:
        case 1:
            hugePages.setMode(value);
            return 0;
        default:
            return 1;
        }
#else
    suppress_unused_warning(param);
    suppress_unused_warning(value);
#endif
    return 1;
}
