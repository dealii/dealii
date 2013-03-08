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
#include <errno.h>

//! Define the main synchronization method
/** It should be specified before including LifoList.h */
#define FINE_GRAIN_LOCKS
#include "LifoList.h"

#if USE_PTHREAD
    #define TlsSetValue_func pthread_setspecific
    #define TlsGetValue_func pthread_getspecific
    typedef pthread_key_t tls_key_t;
    #include <sched.h>
    inline void do_yield() {sched_yield();}

#elif USE_WINTHREAD
    #define _WIN32_WINNT 0x0400
    #include "tbb/machine/windows_api.h"
    #define TlsSetValue_func TlsSetValue
    #define TlsGetValue_func TlsGetValue
    typedef DWORD tls_key_t;
    inline void do_yield() {SwitchToThread();}

#else
    #error Must define USE_PTHREAD or USE_WINTHREAD

#endif


#define FREELIST_NONBLOCKING 1

void mallocThreadShutdownNotification(void* arg);

namespace rml {
namespace internal {

class Block;

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
const unsigned int startupAllocObjSizeMark = ~(unsigned int)0;

/*
 * This number of bins in the TLS that leads to blocks that we can allocate in.
 */
const uint32_t numBlockBinLimit = 31;

/*
 * The following constant is used to define the size of struct Block, the block header.
 * The intent is to have the size of a Block multiple of the cache line size, this allows us to
 * get good alignment at the cost of some overhead equal to the amount of padding included in the Block.
 */

const int blockHeaderAlignment = 64; // a common size of a cache line


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
    void *allocate(size_t size);
    void free(void* ptr);
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
struct TLSData;

class LocalBlockFields : public BlockI {
protected:
    Block       *next;            /* This field needs to be on a 16K boundary and the first field in the block
                                     for LIFO lists to work. */
    Block       *previous;        /* Use double linked list to speed up removal */
    unsigned int objectSize;
    ThreadId     owner;
    FreeObject  *bumpPtr;         /* Bump pointer moves from the end to the beginning of a block */
    FreeObject  *freeList;
    BackRefIdx   backRefIdx;
    unsigned int allocatedCount;  /* Number of objects allocated (obviously by the owning thread) */
    bool         isFull;

    friend void *BootStrapBlocks::allocate(size_t size);
    friend class FreeBlockPool;
    friend class StartupBlock;
    friend void BlockI::initialize(void *bumpPtr);
};
    
class Block : public LocalBlockFields {
    size_t       __pad_local_fields[(blockHeaderAlignment-sizeof(LocalBlockFields))/sizeof(size_t)];
    FreeObject  *publicFreeList;
    Block       *nextPrivatizable;
    size_t       __pad_public_fields[(blockHeaderAlignment-2*sizeof(void*))/sizeof(size_t)];

public:
    static Block *getEmpty(size_t size);

    inline FreeObject* allocate();
    inline FreeObject *allocateFromFreeList();
    inline bool emptyEnoughToUse();
    bool freeListNonNull() { return freeList; }
    void freePublicObject(FreeObject *objectToFree);
    inline void freeOwnObject(FreeObject *objectToFree);
    void returnEmpty(bool poolTheBlock);
    void privatizePublicFreeList();
    void restoreBumpPtr();
    void privatizeOrphaned(Bin *bin);
    void shareOrphaned(const Bin *bin);
    unsigned int getSize() const { return objectSize; }
    const BackRefIdx *getBackRefIdx() const { return &backRefIdx; }
    bool ownBlock() const { return owner.own(); }
    bool isStartupAllocObject() const { return objectSize == startupAllocObjSizeMark; }
    inline FreeObject *findObjectToFree(void *object) const;
    bool checkFreePrecond() const { return allocatedCount>0; }
    const BackRefIdx *getBackRef() const { return &backRefIdx; }
    
protected:
    static Block *getRaw(bool startup);
    void cleanBlockHeader();

private:
    static const float emptyEnoughRatio; /* "Reactivate" a block if this share of its objects is free. */

    inline FreeObject *allocateFromBumpPtr();
    void initEmptyBlock(size_t size);
    inline FreeObject *findAllocatedObject(const void *address) const;
    inline bool isProperlyPlaced(const void *object) const;

    friend class Bin;
    friend void ::mallocThreadShutdownNotification(void* arg);
    friend BlockI *BlockI::getRawBlock(bool startup);
};

const float Block::emptyEnoughRatio = 1.0 / 4.0;

class Bin {
    Block      *activeBlk;
    Block      *mailbox;
    MallocMutex mailLock;

    static TLSData* createTLS();
public:
    static inline Bin* getAllocationBin(size_t size);

    inline Block* getActiveBlock() const { return activeBlk; }
    inline void setActiveBlock(Block *block);
    inline Block* setPreviousBlockActive();
    Block* getPublicFreeListBlock();
    void moveBlockToBinFront(Block *block);
    void processLessUsedBlock(Block *block);

    void outofTLSBin (Block* block);
    void verifyTLSBin (size_t size) const;
    void pushTLSBin(Block* block);

    friend void ::mallocThreadShutdownNotification(void* arg);
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
 * And there are 5 bins with the following allocation sizes: 1792, 2688, 3968, 5376, 8064.
 * They selected to fit 9, 6, 4, 3, and 2 sizes per a block, and also are multiples of 128.
 * If sizeof(Block) changes from 128, these sizes require close attention!
 */
const uint32_t minFittingIndex = minSegregatedObjectIndex+numSegregatedObjectBins;
const uint32_t numFittingBins = 5;

const uint32_t fittingAlignment = 128;

#define SET_FITTING_SIZE(N) ( (blockSize-sizeof(Block))/N ) & ~(fittingAlignment-1)
const uint32_t fittingSize1 = SET_FITTING_SIZE(9);
const uint32_t fittingSize2 = SET_FITTING_SIZE(6);
const uint32_t fittingSize3 = SET_FITTING_SIZE(4);
const uint32_t fittingSize4 = SET_FITTING_SIZE(3);
const uint32_t fittingSize5 = SET_FITTING_SIZE(2);
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
};

static char globalBinSpace[sizeof(LifoList)*numBlockBinLimit];
static OrphanedBlocks *orphanedBlocks = (OrphanedBlocks*)globalBinSpace;

/*
 * Per-thread pool of 16KB blocks. Idea behind it is to not share with other 
 * threads memory that are likely in local cache(s) of our CPU.
 */
class FreeBlockPool {
    Block *head;
    Block *tail;
    int    size;
    void insertBlock(Block *block);
public:
    static const int POOL_HIGH_MARK = 32;
    static const int POOL_LOW_MARK  = 8;

    Block *getBlock();
    void returnBlock(Block *block);
    void releaseAllBlocks();
};

struct TLSData {
    Bin           bin[numBlockBinLimit];
    FreeBlockPool pool;
};

#if MALLOC_CHECK_RECURSION
MallocMutex RecursiveMallocCallProtector::rmc_mutex;
pthread_t   RecursiveMallocCallProtector::owner_thread;
void       *RecursiveMallocCallProtector::autoObjPtr;
bool        RecursiveMallocCallProtector::mallocRecursionDetected;
#if __FreeBSD__
bool        RecursiveMallocCallProtector::canUsePthread;
#endif

#endif

/*********** Code to provide thread ID and a thread-local void pointer **********/

static tls_key_t TLS_pointer_key;

static inline TLSData* getThreadMallocTLS() {
    TLSData *result;
    result = (TLSData *)TlsGetValue_func( TLS_pointer_key );
// The assert below is incorrect: with lazy initialization, it fails on the first call of the function.
//    MALLOC_ASSERT( result, "Memory allocator not initialized" );
    return result;
}

static inline void  setThreadMallocTLS( TLSData * newvalue ) {
    RecursiveMallocCallProtector scoped;
    TlsSetValue_func( TLS_pointer_key, newvalue );
}

/*********** End code to provide thread ID and a TLS pointer **********/

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
 */
#if _WIN64 && _MSC_VER>=1400 && !__INTEL_COMPILER
extern "C" unsigned char _BitScanReverse( unsigned long* i, unsigned long w );
#pragma intrinsic(_BitScanReverse)
#endif
static inline unsigned int highestBitPos(unsigned int n)
{
    unsigned int pos;
#if __ARCH_x86_32||__ARCH_x86_64

# if __linux__||__APPLE__||__FreeBSD__||__sun||__MINGW32__
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

#elif __ARCH_ipf || __ARCH_other
    static unsigned int bsr[16] = {0,6,7,7,8,8,8,8,9,9,9,9,9,9,9,9};
    MALLOC_ASSERT( n>=64 && n<1024, ASSERT_TEXT );
    pos = bsr[ n>>6 ];
#else
#   error highestBitPos() not implemented for this platform
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

/*
 * Initialization code.
 *
 */

/*
 * Forward Refs
 */

static BootStrapBlocks bootStrapBlocks;

void *BootStrapBlocks::allocate(size_t size)
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
                bootStrapBlock = Block::getEmpty(size);
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

/********* Thread and block related code      *************/

void Bin::verifyTLSBin (size_t size) const
{
#if MALLOC_DEBUG
/* The debug version verifies the TLSBin as needed */
    Bin*   tlsBin = getThreadMallocTLS()->bin;
    uint32_t index = getIndex(size);
    uint32_t objSize = getObjectSize(size);

    MALLOC_ASSERT( this == tlsBin+index, ASSERT_TEXT );

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

TLSData* Bin::createTLS()
{
    MALLOC_ASSERT( sizeof(TLSData) >= sizeof(Bin) * numBlockBins + sizeof(FreeBlockPool), ASSERT_TEXT );
    TLSData* tls = (TLSData*) bootStrapBlocks.allocate(sizeof(TLSData));
    if ( !tls ) return NULL;
    /* the block contains zeroes after bootStrapMalloc, so bins are initialized */
#if MALLOC_DEBUG
    for (uint32_t i = 0; i < numBlockBinLimit; i++) {
        MALLOC_ASSERT( tls->bin[i].activeBlk == 0, ASSERT_TEXT );
        MALLOC_ASSERT( tls->bin[i].mailbox == 0, ASSERT_TEXT );
    }
#endif
    setThreadMallocTLS(tls);
    return tls;
}

/*
 * Return the bin for the given size. If the TLS bin structure is absent, create it.
 */
Bin* Bin::getAllocationBin(size_t size)
{
    TLSData* tls = getThreadMallocTLS();
    if( !tls ) 
        tls = createTLS();
    MALLOC_ASSERT( tls, ASSERT_TEXT );
    return tls->bin + getIndex(size);
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

BlockI *BlockI::getRawBlock(bool startup) { return Block::getRaw(startup); }

void BlockI::initialize(void *ptr) { ((LocalBlockFields*)this)->bumpPtr = (FreeObject*)ptr; }

bool Block::emptyEnoughToUse()
{
    const float threshold = (blockSize - sizeof(Block)) * (1-emptyEnoughRatio);

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
    bumpPtr = (FreeObject *)((uintptr_t)this + blockSize - objectSize);
    freeList = NULL;
    isFull = 0;
}

void Block::freeOwnObject(FreeObject *objectToFree)
{
    objectToFree->next = freeList;
    freeList = objectToFree;
    allocatedCount--;
    MALLOC_ASSERT( allocatedCount < (blockSize-sizeof(Block))/objectSize, ASSERT_TEXT );
#if COLLECT_STATISTICS
    if (getActiveBlock(getAllocationBin(block->objectSize)) != block)
        STAT_increment(myTid, getIndex(block->objectSize), freeToInactiveBlock);
    else
        STAT_increment(myTid, getIndex(block->objectSize), freeToActiveBlock);
#endif
    if (isFull) {
        if (emptyEnoughToUse())
            Bin::getAllocationBin(objectSize)->moveBlockToBinFront(this);
    } else {
        if (allocatedCount==0 && publicFreeList==NULL)
            Bin::getAllocationBin(objectSize)->processLessUsedBlock(this);
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
        MALLOC_ASSERT( allocatedCount <= (blockSize-sizeof(Block))/objectSize, ASSERT_TEXT );
        /* other threads did not change the counter freeing our blocks */
        allocatedCount--;
        while( isSolidPtr(temp->next) ){ // the list will end with either NULL or UNUSABLE
            temp = temp->next;
            allocatedCount--;
        }
        MALLOC_ASSERT( allocatedCount < (blockSize-sizeof(Block))/objectSize, ASSERT_TEXT );
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

    publicFreeList = NULL;
}

void Block::initEmptyBlock(size_t size)
{
    // Having getIndex and getObjectSize called next to each other
    // allows better compiler optimization as they basically share the code.
    unsigned int index = getIndex(size);
    unsigned int objSz = getObjectSize(size); 
    Bin* tlsBin = getThreadMallocTLS()->bin;

    cleanBlockHeader();
    objectSize = objSz;
    owner = ThreadId::get();
    // bump pointer should be prepared for first allocation - thus mode it down to objectSize
    bumpPtr = (FreeObject *)((uintptr_t)this + blockSize - objectSize);

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
    result = (Block *) bins[index].pop();
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
    bins[index].push((void **)block);
}

void FreeBlockPool::insertBlock(Block *block)
{
    size++;
    block->next = head;
    head = block;
    if (!tail)
        tail = block;
}

Block *FreeBlockPool::getBlock()
{
    Block *result = head;
    if (head) {
        size--;
        head = head->next;
        if (!head)
            tail = NULL;
    }
    return result;
}

void FreeBlockPool::returnBlock(Block *block)
{
    MALLOC_ASSERT( size <= POOL_HIGH_MARK, ASSERT_TEXT );
    if (size == POOL_HIGH_MARK) {
        // release cold blocks and add hot one
        Block *headToFree = head, 
              *tailToFree = tail;
        for (int i=0; i<POOL_LOW_MARK-2; i++)
            headToFree = headToFree->next;
        tail = headToFree;
        headToFree = headToFree->next;
        tail->next = NULL;
        size = POOL_LOW_MARK-1;
        for (Block *currBl = headToFree; currBl; currBl = currBl->next)
            removeBackRef(currBl->backRefIdx);
        freeBlocks.putList(headToFree, tailToFree);
    }
    insertBlock(block);
}

void FreeBlockPool::releaseAllBlocks()
{
    if (head) {
        for (Block *currBl = head; currBl; currBl = currBl->next)
            removeBackRef(currBl->backRefIdx);
        freeBlocks.putList(head, tail);
    }
}

/* Return an empty uninitialized block in a non-blocking fashion. */
Block *Block::getRaw(bool startup)
{
    Block *result = NULL;
    Block *bigBlock;

    if (! (bigBlock = static_cast<Block*>(freeBlocks.get(startup)))) return NULL;

    // check alignment
    MALLOC_ASSERT( isAligned( bigBlock, blockSize ), ASSERT_TEXT );
    MALLOC_ASSERT( isAligned( bigBlock->bumpPtr, blockSize ), ASSERT_TEXT );
    // block should be at least as big as blockSize; otherwise the previous block can be damaged.
    MALLOC_ASSERT( (uintptr_t)bigBlock->bumpPtr >= (uintptr_t)bigBlock + blockSize, ASSERT_TEXT );
    bigBlock->bumpPtr = (FreeObject *)((uintptr_t)bigBlock->bumpPtr - blockSize);
    result = (Block *)bigBlock->bumpPtr;
    if ( result!=bigBlock ) {
        TRACEF(( "[ScalableMalloc trace] Pushing partial rest of block back on.\n" ));
        freeBlocks.put(bigBlock, startup);
    }

    return result;
}

/* Return an empty uninitialized block in a non-blocking fashion. */
Block *Block::getEmpty(size_t size)
{
    Block *result = NULL;
    TLSData* tls = getThreadMallocTLS();
    if (tls)
        result = tls->pool.getBlock();
    if (!result) {
        BackRefIdx backRefIdx = BackRefIdx::newBackRef(/*largeObj=*/false);
        if (backRefIdx.isInvalid() || !(result = getRaw(/*startup=*/false)))
            return NULL;
        setBackRef(backRefIdx, result);
        result->backRefIdx = backRefIdx;
    }
    if (result) {
        result->initEmptyBlock(size);
        STAT_increment(result->owner, getIndex(result->objectSize), allocBlockNew);
    }
    return result;
}

/* We have a block give it back to the malloc block manager */
void Block::returnEmpty(bool poolTheBlock)
{
    // it is caller's responsibility to ensure no data is lost before calling this
    MALLOC_ASSERT( allocatedCount==0, ASSERT_TEXT );
    MALLOC_ASSERT( publicFreeList==NULL, ASSERT_TEXT );
    MALLOC_ASSERT( !poolTheBlock || next == NULL, ASSERT_TEXT );
    MALLOC_ASSERT( !poolTheBlock || previous == NULL, ASSERT_TEXT );
    STAT_increment(owner, getIndex(objectSize), freeBlockBack);

    cleanBlockHeader();

    nextPrivatizable = NULL;

    objectSize = 0;
    owner.invalid();
    // for an empty block, bump pointer should point right after the end of the block
    bumpPtr = (FreeObject *)((uintptr_t)this + blockSize);
    if (poolTheBlock) {
        MALLOC_ASSERT(getThreadMallocTLS(), "Is TLS still not initialized?");
        getThreadMallocTLS()->pool.returnBlock(this);
    }
    else {
        removeBackRef(backRefIdx);
        freeBlocks.put(this, /*startup=*/false);
    }
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

#if MALLOC_CHECK_RECURSION

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
    size_t availableSize() {
        return blockSize - ((uintptr_t)bumpPtr - (uintptr_t)this);
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

    StartupBlock *block = (StartupBlock *)getRaw(/*startup=*/true);
    if (!block) return NULL;

    block->cleanBlockHeader();
    setBackRef(backRefIdx, block);
    block->backRefIdx = backRefIdx;
    // use startupAllocObjSizeMark to mark objects from startup block marker
    block->objectSize = startupAllocObjSizeMark;
    block->bumpPtr = (FreeObject *)((uintptr_t)block + sizeof(StartupBlock));
    return block;
}

/* TODO: Function is called when malloc nested call is detected, so simultaneous
   usage from different threads are unprobable, so block pre-allocation 
   can be not useful, and the code might be simplified. */
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
    // speculatively allocates newBlock to later use or return it as unused
    if (!firstStartupBlock || firstStartupBlock->availableSize() < reqSize)
        if (!(newBlock = StartupBlock::getBlock()))
            return NULL;

    {
        MallocMutex::scoped_lock scoped_cs(startupMallocLock);
    
        if (!firstStartupBlock || firstStartupBlock->availableSize() < reqSize) {
            if (!newBlock && !(newBlock = StartupBlock::getBlock()))
                return NULL;
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
        newBlock->returnEmpty(/*poolTheBlock=*/false);

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
                      && (uintptr_t)ptr+StartupBlock::msize(ptr)<=(uintptr_t)this+blockSize, 
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
        blockToRelease->returnEmpty(/*poolTheBlock=*/false);
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

// TODO: add error handling, and on error do something better than exit(1)
    if (!initBackRefMaster() || !freeBlocks.bootstrap(NULL, NULL, 0)) {
        fprintf (stderr, "The memory manager cannot access sufficient memory to initialize; exiting \n");
        exit(1);
    }
// Create keys for thread-local storage and for thread id
#if USE_WINTHREAD
    TLS_pointer_key = TlsAlloc();
#else
    int status1 = pthread_key_create( &TLS_pointer_key, mallocThreadShutdownNotification );
    if ( status1 ) {
        fprintf (stderr, "The memory manager cannot create tls key during initialization; exiting \n");
        exit(1);
    }
#endif /* USE_WINTHREAD */
    ThreadId::init();
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
    MALLOC_ASSERT( allocatedCount < (blockSize-sizeof(Block))/objectSize, ASSERT_TEXT );
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
        MALLOC_ASSERT( allocatedCount < (blockSize-sizeof(Block))/objectSize, ASSERT_TEXT );
        allocatedCount++;
        STAT_increment(owner, getIndex(objectSize), allocBumpPtrUsed);
    }
    return result;
}

inline FreeObject* Block::allocate()
{
    FreeObject *result;

    MALLOC_ASSERT( owner.own(), ASSERT_TEXT );

    /* for better cache locality, first looking in the free list. */
    if ( (result = allocateFromFreeList()) ) {
        return result;
    }
    MALLOC_ASSERT( !freeList, ASSERT_TEXT );

    /* if free list is empty, try thread local bump pointer allocation. */
    if ( (result = allocateFromBumpPtr()) ) {
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

void Bin::processLessUsedBlock(Block *block)
{
    if (block != activeBlk) {
        /* We are not actively using this block; return it to the general block pool */
        outofTLSBin(block);
        block->returnEmpty(/*poolTheBlock=*/true);
    } else {
        /* all objects are free - let's restore the bump pointer */
        block->restoreBumpPtr();
    }
}

/*
 * All aligned allocations fall into one of the following categories:
 *  1. if both request size and alignment are <= maxSegregatedObjectSize,
 *       we just align the size up, and request this amount, because for every size
 *       aligned to some power of 2, the allocated object is at least that aligned.
 * 2. for bigger size, check if already guaranteed fittingAlignment is enough.
 * 3. if size+alignment<minLargeObjectSize, we take an object of fittingSizeN and align
 *       its address up; given such pointer, scalable_free could find the real object.
 * 4. otherwise, aligned large object is allocated.
 */
static void *allocateAligned(size_t size, size_t alignment)
{
    MALLOC_ASSERT( isPowerOfTwo(alignment), ASSERT_TEXT );

    void *result;
    if (size<=maxSegregatedObjectSize && alignment<=maxSegregatedObjectSize)
        result = scalable_malloc(alignUp(size? size: sizeof(size_t), alignment));
    else if (size<minLargeObjectSize && alignment<=fittingAlignment)
        result = scalable_malloc(size);
    else if (size+alignment < minLargeObjectSize) {
        void *unaligned = scalable_malloc(size+alignment);
        if (!unaligned) return NULL;
        result = alignUp(unaligned, alignment);
    } else {
        /* This can be the first allocation call. */
        if (!isMallocInitialized()) 
            doInitialization();
        // take into account only alignment that are higher then natural
        result = mallocLargeObject(size, largeObjectAlignment>alignment? 
                                         largeObjectAlignment: alignment);
    }

    MALLOC_ASSERT( isAligned(result, alignment), ASSERT_TEXT );
    return result;
}

static void *reallocAligned(void *ptr, size_t size, size_t alignment = 0)
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
            result = alignment ? allocateAligned(size, alignment) : scalable_malloc(size);
        }
    } else {
        Block* block = (Block *)alignDown(ptr, blockSize);
        copySize = block->getSize();
        if (size <= copySize && (0==alignment || isAligned(ptr, alignment))) {
            return ptr;
        } else {
            result = alignment ? allocateAligned(size, alignment) : scalable_malloc(size);
        }
    }
    if (result) {
        memcpy(result, ptr, copySize<size? copySize: size);
        scalable_free(ptr);
    }
    return result;
}

/* A predicate checks if an object is properly placed inside its block */
inline bool Block::isProperlyPlaced(const void *object) const
{
    return 0 == ((uintptr_t)this + blockSize - (uintptr_t)object) % objectSize;
}

/* Finds the real object inside the block */
FreeObject *Block::findAllocatedObject(const void *address) const
{
    // calculate offset from the end of the block space
    uintptr_t offset = (uintptr_t)this + blockSize - (uintptr_t)address;
    MALLOC_ASSERT( offset<blockSize-sizeof(Block), ASSERT_TEXT );
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

    return idx.isLargeObject() && getBackRef(idx) == header;
}

static inline bool isSmallObject (void *ptr)
{
    void* expected = alignDown(ptr, blockSize);
    const BackRefIdx* idx = ((Block*)expected)->getBackRef();

    return expected == getBackRef(safer_dereference(idx));
}

/**** Check if an object was allocated by scalable_malloc ****/
static inline bool isRecognized (void* ptr)
{
    return isLargeObject(ptr) || isSmallObject(ptr);
}

static inline void freeSmallObject (void *object) 
{
    /* mask low bits to get the block */
    Block *block = (Block *)alignDown(object, blockSize);
    MALLOC_ASSERT( block->checkFreePrecond(), ASSERT_TEXT );

#if MALLOC_CHECK_RECURSION
    if (block->isStartupAllocObject()) {
        ((StartupBlock *)block)->free(object);
        return;
    }
#endif
    FreeObject *objectToFree = block->findObjectToFree(object);

    if (block->ownBlock())
        block->freeOwnObject(objectToFree);
    else /* Slower path to add to the shared list, the allocatedCount is updated by the owner thread in malloc. */
        block->freePublicObject(objectToFree);

}

} // namespace internal
} // namespace rml

using namespace rml::internal;

/*
 * When a thread is shutting down this routine should be called to remove all the thread ids
 * from the malloc blocks and replace them with a NULL thread id.
 *
 */
#if MALLOC_TRACE
static unsigned int threadGoingDownCount = 0;
#endif

/*
 * for pthreads, the function is set as a callback in pthread_key_create for TLS bin.
 * it will be automatically called at thread exit with the key value as the argument.
 *
 * for Windows, it should be called directly e.g. from DllMain; the argument can be NULL
 * one should include "TypeDefinitions.h" for the declaration of this function.
*/
extern "C" void mallocThreadShutdownNotification(void* arg)
{
    TLSData *tls;
    Block *threadBlock;
    Block *threadlessBlock;
    unsigned int index;

    // Check whether TLS has been initialized
    if (!isMallocInitialized()) return;

    TRACEF(( "[ScalableMalloc trace] Thread id %d blocks return start %d\n",
             getThreadId(),  threadGoingDownCount++ ));
#ifdef USE_WINTHREAD
    tls = getThreadMallocTLS();
#else
    tls = (TLSData*)arg;
#endif
    if (tls) {
        Bin *tlsBin = tls->bin;
        tls->pool.releaseAllBlocks();

        for (index = 0; index < numBlockBins; index++) {
            if (tlsBin[index].activeBlk==NULL)
                continue;
            threadlessBlock = tlsBin[index].activeBlk->previous;
            while (threadlessBlock) {
                threadBlock = threadlessBlock->previous;
                if (threadlessBlock->allocatedCount==0 && threadlessBlock->publicFreeList==NULL) {
                    /* we destroy the thread, so not use its block pool */
                    threadlessBlock->returnEmpty(/*poolTheBlock=*/false);
                } else {
                    orphanedBlocks->put(tlsBin+index, threadlessBlock);
                }
                threadlessBlock = threadBlock;
            }
            threadlessBlock = tlsBin[index].activeBlk;
            while (threadlessBlock) {
                threadBlock = threadlessBlock->next;
                if (threadlessBlock->allocatedCount==0 && threadlessBlock->publicFreeList==NULL) {
                    /* we destroy the thread, so not use its block pool */
                    threadlessBlock->returnEmpty(/*poolTheBlock=*/false);
                } else {
                    orphanedBlocks->put(tlsBin+index, threadlessBlock);
                }
                threadlessBlock = threadBlock;
            }
            tlsBin[index].activeBlk = 0;
        }
        bootStrapBlocks.free(tls);
        setThreadMallocTLS(NULL);
    }

    TRACEF(( "[ScalableMalloc trace] Thread id %d blocks return end\n", getThreadId() ));
}

extern "C" void mallocProcessShutdownNotification(void)
{
#if COLLECT_STATISTICS
    ThreadId nThreads = ThreadIdCount;
    for( int i=1; i<=nThreads && i<MAX_THREADS; ++i )
        STAT_print(i);
#endif
}

/********* The malloc code          *************/

extern "C" void * scalable_malloc(size_t size)
{
    Bin* bin;
    Block * mallocBlock;
    FreeObject *result = NULL;

    if (!size) size = sizeof(size_t);

#if MALLOC_CHECK_RECURSION
    if (RecursiveMallocCallProtector::sameThreadActive()) {
        result = size<minLargeObjectSize? StartupBlock::allocate(size) : 
              (FreeObject*)mallocLargeObject(size, blockSize, /*startupAlloc=*/ true);
        if (!result) errno = ENOMEM;
        return result;
    }
#endif

    if (!isMallocInitialized()) 
        doInitialization();

    /*
     * Use Large Object Allocation
     */
    if (size >= minLargeObjectSize) {
        result = (FreeObject*)mallocLargeObject(size, largeObjectAlignment);
        if (!result) errno = ENOMEM;
        return result;
    }

    /*
     * Get an element in thread-local array corresponding to the given size;
     * It keeps ptr to the active block for allocations of this size
     */
    bin = Bin::getAllocationBin(size);
    if ( !bin ) {
        errno = ENOMEM;
        return NULL;
    }

    /* Get the block of you want to try to allocate in. */
    mallocBlock = bin->getActiveBlock();

    if (mallocBlock) {
        do {
            if( (result = mallocBlock->allocate()) ) {
                return result;
            }
            // the previous block, if any, should be empty enough
        } while( (mallocBlock = bin->setPreviousBlockActive()) );
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
        if ( (result = mallocBlock->allocateFromFreeList()) ) {
            return result;
        }
        /* Else something strange happened, need to retry from the beginning; */
        TRACEF(( "[ScalableMalloc trace] Something is wrong: no objects in public free list; reentering.\n" ));
        return scalable_malloc(size);
    }

    /*
     * no suitable own blocks, try to get a partial block that some other thread has discarded.
     */
    mallocBlock = orphanedBlocks->get(bin, size);
    while (mallocBlock) {
        bin->pushTLSBin(mallocBlock);
        bin->setActiveBlock(mallocBlock); // TODO: move under the below condition?
        if( (result = mallocBlock->allocate()) ) {
            return result;
        }
        mallocBlock = orphanedBlocks->get(bin, size);
    }

    /*
     * else try to get a new empty block
     */
    mallocBlock = Block::getEmpty(size);
    if (mallocBlock) {
        bin->pushTLSBin(mallocBlock);
        bin->setActiveBlock(mallocBlock);
        if( (result = mallocBlock->allocate()) ) {
            return result;
        }
        /* Else something strange happened, need to retry from the beginning; */
        TRACEF(( "[ScalableMalloc trace] Something is wrong: no objects in empty block; reentering.\n" ));
        return scalable_malloc(size);
    }
    /*
     * else nothing works so return NULL
     */
    TRACEF(( "[ScalableMalloc trace] No memory found, returning NULL.\n" ));
    errno = ENOMEM;
    return NULL;
}

/********* End the malloc code      *************/

/********* The free code            *************/

extern "C" void scalable_free (void *object) {
if (!object)
        return;

    MALLOC_ASSERT(isRecognized(object), "Invalid pointer in scalable_free detected.");

    if (isLargeObject(object))
        freeLargeObject(object);
    else
        freeSmallObject(object);
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
    if (isLargeObject(object))
        freeLargeObject(object);
    else if (isSmallObject(object))
        freeSmallObject(object);
    else if (original_free)
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
    /* corner cases left out of reallocAligned to not deal with errno there */
    if (!ptr) {
        return scalable_malloc(size);
    }
    if (!size) {
        scalable_free(ptr);
        return NULL;
    }
    void* tmp = reallocAligned(ptr, size, 0);
    if (!tmp) errno = ENOMEM;
    return tmp;
}

/*
 * A variant that provides additional memory safety, by checking whether the given address
 * was obtained with this allocator, and if not redirecting to the provided alternative call.
 */
extern "C" void* safer_scalable_realloc (void* ptr, size_t sz, void* original_realloc) 
{
    if (!ptr) {
        return scalable_malloc(sz);
    }
    if (isRecognized(ptr)) {
        if (!sz) {
            scalable_free(ptr);
            return NULL;
        }
        void* tmp = reallocAligned(ptr, sz, 0);
        if (!tmp) errno = ENOMEM;
        return tmp;
    }
#if USE_WINTHREAD
    else if (original_realloc && sz) {
            orig_ptrs *original_ptrs = static_cast<orig_ptrs*>(original_realloc);
            if ( original_ptrs->orig_msize ){
                size_t oldSize = original_ptrs->orig_msize(ptr);
                void *newBuf = scalable_malloc(sz);
                if (newBuf) {
                    memcpy(newBuf, ptr, sz<oldSize? sz : oldSize);
                    if ( original_ptrs->orig_free ){
                        original_ptrs->orig_free( ptr );
                    }
                }
                return newBuf;
             }
    }
#else
    else if (original_realloc) {
        typedef void* (*realloc_ptr_t)(void*,size_t);
        realloc_ptr_t original_realloc_ptr;
        (void *&)original_realloc_ptr = original_realloc;
        return original_realloc_ptr(ptr,sz);
    }
#endif
    return NULL;
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
    void* result = scalable_malloc(arraySize);
    if (result)
        memset(result, 0, arraySize);
    return result;
}

/********* End code for scalable_calloc   ***********/

/********* Code for aligned allocation API **********/

extern "C" int scalable_posix_memalign(void **memptr, size_t alignment, size_t size)
{
    if ( !isPowerOfTwoMultiple(alignment, sizeof(void*)) )
        return EINVAL;
    void *result = allocateAligned(size, alignment);
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
    void* tmp = allocateAligned(size, alignment);
    if (!tmp) 
        errno = ENOMEM;
    return tmp;
}

extern "C" void * scalable_aligned_realloc(void *ptr, size_t size, size_t alignment)
{
    /* corner cases left out of reallocAligned to not deal with errno there */
    if (!isPowerOfTwo(alignment)) {
        errno = EINVAL;
        return NULL;
    }
    if (!ptr) {
        return allocateAligned(size, alignment);
    }
    if (!size) {
        scalable_free(ptr);
        return NULL;
    }

    void* tmp = reallocAligned(ptr, size, alignment);
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
    if (!ptr) {
        return allocateAligned(size, alignment);
    }
    if (isRecognized(ptr)) {
        if (!size) {
            scalable_free(ptr);
            return NULL;
        }
        void* tmp = reallocAligned(ptr, size, alignment);
        if (!tmp) errno = ENOMEM;
        return tmp;
    }
#if USE_WINTHREAD
    else {
        orig_ptrs *original_ptrs = static_cast<orig_ptrs*>(orig_function);
        if (size) {
            if ( original_ptrs->orig_msize ){
                size_t oldSize = original_ptrs->orig_msize(ptr);
                void *newBuf = allocateAligned(size, alignment);
                if (newBuf) {
                    memcpy(newBuf, ptr, size<oldSize? size : oldSize);
                    if ( original_ptrs->orig_free ){
                        original_ptrs->orig_free( ptr );
                    }
                }
                return newBuf;
            }else{
                //We can't do anything with this. Just keeping old pointer
                return NULL;
            }
        } else {
            if ( original_ptrs->orig_free ){
                original_ptrs->orig_free( ptr );
            }
            return NULL;
        }
    }
#endif
    return NULL;
}

extern "C" void scalable_aligned_free(void *ptr)
{
    scalable_free(ptr);
}

/********* end code for aligned allocation API **********/

/********* Code for scalable_msize       ***********/

/*
 * Returns the size of a memory block allocated in the heap.
 */
extern "C" size_t scalable_msize(void* ptr)
{
    if (ptr) {
        MALLOC_ASSERT(isRecognized(ptr), "Invalid pointer in scalable_msize detected.");
        if (isLargeObject(ptr)) {
            LargeMemoryBlock* lmb = ((LargeObjectHdr*)ptr - 1)->memoryBlock;
            return lmb->objectSize;
        } else {
            Block* block = (Block *)alignDown(ptr, blockSize);
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

/*
 * A variant that provides additional memory safety, by checking whether the given address
 * was obtained with this allocator, and if not redirecting to the provided alternative call.
 */
extern "C" size_t safer_scalable_msize (void *object, size_t (*original_msize)(void*)) 
{
    if (object) {
        // Check if the memory was allocated by scalable_malloc
        if (isRecognized(object))
            return scalable_msize(object);
        else if (original_msize)
            return original_msize(object);
    }
    // object is NULL or unknown
    errno = EINVAL;
    return 0;
}

/********* End code for scalable_msize   ***********/
