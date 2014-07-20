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

#include <string.h>   /* for memset */
#include <errno.h>
#include "tbbmalloc_internal.h"

namespace rml {
namespace internal {

/*********** Code to acquire memory from the OS or other executive ****************/

/*
  syscall/malloc can set non-zero errno in case of failure,
  but later allocator might be able to find memory to fulfil the request.
  And we do not want changing of errno by successful scalable_malloc call.
  To support this, restore old errno in (get|free)RawMemory, and set errno
  in frontend just before returning to user code.
  Please note: every syscall/libc call used inside scalable_malloc that
  sets errno must be protected this way, not just memory allocation per se.
*/

#if USE_DEFAULT_MEMORY_MAPPING
#include "MapMemory.h"
#else
/* assume MapMemory and UnmapMemory are customized */
#endif

void* getRawMemory (size_t size, bool hugePages) {
    return MapMemory(size, hugePages);
}

bool freeRawMemory (void *object, size_t size) {
    return UnmapMemory(object, size);
}

void HugePagesStatus::registerAllocation(bool gotPage)
{
    if (gotPage) {
        if (!wasObserved)
            FencedStore(wasObserved, 1);
    } else
        FencedStore(enabled, 0);
    // reports huge page status only once
    if (needActualStatusPrint
        && AtomicCompareExchange(needActualStatusPrint, 0, 1))
        doPrintStatus(gotPage, "available");
}

void HugePagesStatus::registerReleasing(size_t size)
{
    // We: 1) got huge page at least once,
    // 2) something that looks like a huge page is been released,
    // and 3) user requested huge pages,
    // so a huge page might be available at next allocation.
    // TODO: keep page status in regions and use exact check here
    // Use isPowerOfTwoMultiple because it's faster then generic reminder.
    if (FencedLoad(wasObserved) && isPowerOfTwoMultiple(size, pageSize))
        FencedStore(enabled, requestedMode.get());
}

void HugePagesStatus::printStatus() {
    doPrintStatus(requestedMode.get(), "requested");
    if (requestedMode.get()) { // report actual status iff requested
        if (pageSize)
            FencedStore(needActualStatusPrint, 1);
        else
            doPrintStatus(/*state=*/false, "available");
    }
}

void HugePagesStatus::doPrintStatus(bool state, const char *stateName)
{
    fprintf(stderr, "TBBmalloc: huge pages\t%s%s\n",
            state? "" : "not ", stateName);
}

void *Backend::getRawMem(size_t &size) const
{
    if (extMemPool->userPool()) {
        size = alignUpGeneric(size, extMemPool->granularity);
        return (*extMemPool->rawAlloc)(extMemPool->poolId, size);
    }
    // try to get them at 1st allocation and still use, if successful
    // if 1st try is unsuccessful, no more trying
    if (FencedLoad(hugePages.enabled)) {
        size_t hugeSize = alignUpGeneric(size, hugePages.getSize());
        void *res = getRawMemory(hugeSize, /*hugePages=*/true);
        hugePages.registerAllocation(res);
        if (res) {
            size = hugeSize;
            return res;
        }
    }
    size_t granSize = alignUpGeneric(size, extMemPool->granularity);
    if (void *res = getRawMemory(granSize, /*hugePages=*/false)) {
        size = granSize;
        return res;
    }
    return NULL;
}

void Backend::freeRawMem(void *object, size_t size) const
{
    if (extMemPool->userPool())
        (*extMemPool->rawFree)(extMemPool->poolId, object, size);
    else {
        hugePages.registerReleasing(size);
        freeRawMemory(object, size);
    }
}

/********* End memory acquisition code ********************************/

// Protected object size. After successful locking returns size of locked block,
// and releasing requires setting block size.
class GuardedSize : tbb::internal::no_copy {
    uintptr_t value;
public:
    enum State {
        LOCKED,
        COAL_BLOCK,        // block is coalescing now
        MAX_LOCKED_VAL = COAL_BLOCK,
        LAST_REGION_BLOCK, // used to mark last block in region
        // values after this are "normal" block sizes
        MAX_SPEC_VAL = LAST_REGION_BLOCK
    };

    void initLocked() { value = LOCKED; }
    void makeCoalscing() {
        MALLOC_ASSERT(value == LOCKED, ASSERT_TEXT);
        value = COAL_BLOCK;
    }
    size_t tryLock(State state) {
        size_t szVal, sz;
        MALLOC_ASSERT(state <= MAX_LOCKED_VAL, ASSERT_TEXT);
        for (;;) {
            sz = FencedLoad((intptr_t&)value);
            if (sz <= MAX_LOCKED_VAL)
                break;
            szVal = AtomicCompareExchange((intptr_t&)value, state, sz);

            if (szVal==sz)
                break;
        }
        return sz;
    }
    void unlock(size_t size) {
        MALLOC_ASSERT(value <= MAX_LOCKED_VAL, "The lock is not locked");
        MALLOC_ASSERT(size > MAX_LOCKED_VAL, ASSERT_TEXT);
        FencedStore((intptr_t&)value, size);
    }
    friend void Backend::IndexedBins::verify();
};

struct MemRegion {
    MemRegion *next,    // keep all regions in any pool to release all them on
              *prev;    // pool destroying, 2-linked list to release individual
                        // regions.
    size_t     allocSz, // got from poll callback
               blockSz; // initial and maximal inner block size
    bool       exact;   // region tageted to exact large object allocation
};

// this data must be unmodified while block is in use, so separate it
class BlockMutexes {
protected:
    GuardedSize myL,   // lock for me
                leftL; // lock for left neighbor
};

class FreeBlock : BlockMutexes {
public:
    static const size_t minBlockSize;
    friend void Backend::IndexedBins::verify();

    FreeBlock    *prev,       // in 2-linked list related to bin
                 *next,
                 *nextToFree; // used to form a queue during coalescing
    // valid only when block is in processing, i.e. one is not free and not
    size_t        sizeTmp;    // used outside of backend
    int           myBin;      // bin that is owner of the block
    bool          aligned;
    bool          blockInBin; // this block in myBin already

    FreeBlock *rightNeig(size_t sz) const {
        MALLOC_ASSERT(sz, ASSERT_TEXT);
        return (FreeBlock*)((uintptr_t)this+sz);
    }
    FreeBlock *leftNeig(size_t sz) const {
        MALLOC_ASSERT(sz, ASSERT_TEXT);
        return (FreeBlock*)((uintptr_t)this - sz);
    }

    void initHeader() { myL.initLocked(); leftL.initLocked(); }
    void setMeFree(size_t size) { myL.unlock(size); }
    size_t trySetMeUsed(GuardedSize::State s) { return myL.tryLock(s); }

    void setLeftFree(size_t sz) { leftL.unlock(sz); }
    size_t trySetLeftUsed(GuardedSize::State s) { return leftL.tryLock(s); }

    size_t tryLockBlock() {
        size_t rSz, sz = trySetMeUsed(GuardedSize::LOCKED);

        if (sz <= GuardedSize::MAX_LOCKED_VAL)
            return false;
        rSz = rightNeig(sz)->trySetLeftUsed(GuardedSize::LOCKED);
        if (rSz <= GuardedSize::MAX_LOCKED_VAL) {
            setMeFree(sz);
            return false;
        }
        MALLOC_ASSERT(rSz == sz, ASSERT_TEXT);
        return sz;
    }
    void markCoalescing(size_t blockSz) {
        myL.makeCoalscing();
        rightNeig(blockSz)->leftL.makeCoalscing();
        sizeTmp = blockSz;
        nextToFree = NULL;
    }
    void markUsed() {
        myL.initLocked();
        rightNeig(sizeTmp)->leftL.initLocked();
        nextToFree = NULL;
    }
    static void markBlocks(FreeBlock *fBlock, int num, size_t size) {
        for (int i=1; i<num; i++) {
            fBlock = (FreeBlock*)((uintptr_t)fBlock + size);
            fBlock->initHeader();
        }
    }
};

// Last block in any region. Its "size" field is GuardedSize::LAST_REGION_BLOCK,
// This kind of blocks used to find region header
// and have a possibility to return region back to OS
struct LastFreeBlock : public FreeBlock {
    MemRegion *memRegion;
};

const size_t FreeBlock::minBlockSize = sizeof(FreeBlock);

void CoalRequestQ::putBlock(FreeBlock *fBlock)
{
    MALLOC_ASSERT(fBlock->sizeTmp >= FreeBlock::minBlockSize, ASSERT_TEXT);
    fBlock->markUsed();

    for (;;) {
        FreeBlock *myBlToFree = (FreeBlock*)FencedLoad((intptr_t&)blocksToFree);

        fBlock->nextToFree = myBlToFree;
        if (myBlToFree ==
            (FreeBlock*)AtomicCompareExchange((intptr_t&)blocksToFree,
                                              (intptr_t)fBlock,
                                              (intptr_t)myBlToFree))
            return;
    }
}

FreeBlock *CoalRequestQ::getAll()
{
    for (;;) {
        FreeBlock *myBlToFree = (FreeBlock*)FencedLoad((intptr_t&)blocksToFree);

        if (!myBlToFree)
            return NULL;
        else {
            if (myBlToFree ==
                (FreeBlock*)AtomicCompareExchange((intptr_t&)blocksToFree,
                                                  0, (intptr_t)myBlToFree))
                return myBlToFree;
            else
                continue;
        }
    }
}

// Try to get a block from a bin.
// If the remaining free space would stay in the same bin,
//     split the block without removing it.
// If the free space should go to other bin(s), remove the block.
// alignedBin is true, if all blocks in the bin has slab-aligned right side.
FreeBlock *Backend::IndexedBins::getBlock(int binIdx, BackendSync *sync,
                size_t size, bool needAlignedRes, bool alignedBin, bool wait,
                int *binLocked)
{
    Bin *b = &freeBins[binIdx];
try_next:
    FreeBlock *fBlock = NULL;
    if (b->head) {
        bool locked;
        MallocMutex::scoped_lock scopedLock(b->tLock, wait, &locked);

        if (!locked) {
            if (binLocked) (*binLocked)++;
            return NULL;
        }

        for (FreeBlock *curr = b->head; curr; curr = curr->next) {
            size_t szBlock = curr->tryLockBlock();
            if (!szBlock) {
                goto try_next;
            }

            if (alignedBin || !needAlignedRes) {
                size_t splitSz = szBlock - size;
                // If we got a block as split result,
                // it must have a room for control structures.
                if (szBlock >= size && (splitSz >= FreeBlock::minBlockSize ||
                                        !splitSz))
                    fBlock = curr;
            } else {
                void *newB = alignUp(curr, slabSize);
                uintptr_t rightNew = (uintptr_t)newB + size;
                uintptr_t rightCurr = (uintptr_t)curr + szBlock;
                // appropriate size, and left and right split results
                // are either big enough or non-exitent
                if (rightNew <= rightCurr
                    && (newB==curr ||
                        (uintptr_t)newB-(uintptr_t)curr >= FreeBlock::minBlockSize)
                    && (rightNew==rightCurr ||
                        rightCurr - rightNew >= FreeBlock::minBlockSize))
                    fBlock = curr;
            }
            if (fBlock) {
                // consume must be called before result of removing from a bin
                // is visible externally.
                sync->consume();
                if (alignedBin && needAlignedRes &&
                    Backend::sizeToBin(szBlock-size) == Backend::sizeToBin(szBlock)) {
                    // free remainder of fBlock stay in same bin,
                    // so no need to remove it from the bin
                    // TODO: add more "still here" cases
                    FreeBlock *newFBlock = fBlock;
                    // return block from right side of fBlock
                    fBlock = (FreeBlock*)((uintptr_t)newFBlock + szBlock - size);
                    MALLOC_ASSERT(isAligned(fBlock, slabSize), "Invalid free block");
                    fBlock->initHeader();
                    fBlock->setLeftFree(szBlock - size);
                    newFBlock->setMeFree(szBlock - size);

                    fBlock->sizeTmp = size;
                } else {
                    b->removeBlock(fBlock);
                    if (freeBins[binIdx].empty())
                        bitMask.set(binIdx, false);
                    fBlock->sizeTmp = szBlock;
                }
                break;
            } else { // block size is not valid, search for next block in the bin
                curr->setMeFree(szBlock);
                curr->rightNeig(szBlock)->setLeftFree(szBlock);
            }
        }
    }
    return fBlock;
}

void Backend::Bin::removeBlock(FreeBlock *fBlock)
{
    if (head == fBlock)
        head = fBlock->next;
    if (tail == fBlock)
        tail = fBlock->prev;
    if (fBlock->prev)
        fBlock->prev->next = fBlock->next;
    if (fBlock->next)
        fBlock->next->prev = fBlock->prev;
}

void Backend::IndexedBins::addBlock(int binIdx, FreeBlock *fBlock, size_t blockSz, bool addToTail)
{
    Bin *b = &freeBins[binIdx];

    fBlock->myBin = binIdx;
    fBlock->aligned = toAlignedBin(fBlock, blockSz);
    fBlock->next = fBlock->prev = NULL;
    {
        MallocMutex::scoped_lock scopedLock(b->tLock);
        if (addToTail) {
            fBlock->prev = b->tail;
            b->tail = fBlock;
            if (fBlock->prev)
                fBlock->prev->next = fBlock;
            if (!b->head)
                b->head = fBlock;
        } else {
            fBlock->next = b->head;
            b->head = fBlock;
            if (fBlock->next)
                fBlock->next->prev = fBlock;
            if (!b->tail)
                b->tail = fBlock;
        }
    }
    bitMask.set(binIdx, true);
}

bool Backend::IndexedBins::tryAddBlock(int binIdx, FreeBlock *fBlock, bool addToTail)
{
    bool locked;
    Bin *b = &freeBins[binIdx];

    fBlock->myBin = binIdx;
    fBlock->aligned = toAlignedBin(fBlock, fBlock->sizeTmp);
    if (addToTail) {
        fBlock->next = NULL;
        {
            MallocMutex::scoped_lock scopedLock(b->tLock, /*wait=*/false, &locked);
            if (!locked)
                return false;
            fBlock->prev = b->tail;
            b->tail = fBlock;
            if (fBlock->prev)
                fBlock->prev->next = fBlock;
            if (!b->head)
                b->head = fBlock;
        }
    } else {
        fBlock->prev = NULL;
        {
            MallocMutex::scoped_lock scopedLock(b->tLock, /*wait=*/false, &locked);
            if (!locked)
                return false;
            fBlock->next = b->head;
            b->head = fBlock;
            if (fBlock->next)
                fBlock->next->prev = fBlock;
            if (!b->tail)
                b->tail = fBlock;
        }
    }
    bitMask.set(binIdx, true);
    return true;
}

void Backend::IndexedBins::reset()
{
    for (int i=0; i<Backend::freeBinsNum; i++)
        freeBins[i].reset();
    bitMask.reset();
}

void Backend::IndexedBins::lockRemoveBlock(int binIdx, FreeBlock *fBlock)
{
    MallocMutex::scoped_lock scopedLock(freeBins[binIdx].tLock);
    freeBins[binIdx].removeBlock(fBlock);
    if (freeBins[binIdx].empty())
        bitMask.set(binIdx, false);
}

bool ExtMemoryPool::regionsAreReleaseable() const
{
    return !keepAllMemory && !delayRegsReleasing;
}

// try to allocate num blocks of size Bytes from particular "generic" bin
// needAlignedRes is true if result must be slab-aligned
FreeBlock *Backend::getFromBin(int binIdx, int num, size_t size, bool needAlignedRes,
                               int *binLocked)
{
    FreeBlock *fBlock =
        freeLargeBins.getBlock(binIdx, &bkndSync, num*size, needAlignedRes,
                               /*alignedBin=*/false, /*wait=*/false, binLocked);
    if (fBlock) {
        if (needAlignedRes) {
            size_t fBlockSz = fBlock->sizeTmp;
            uintptr_t fBlockEnd = (uintptr_t)fBlock + fBlockSz;
            FreeBlock *newB = alignUp(fBlock, slabSize);
            FreeBlock *rightPart = (FreeBlock*)((uintptr_t)newB + num*size);

            // Space to use is in the middle,
            // ... return free right part
            if ((uintptr_t)rightPart != fBlockEnd) {
                rightPart->initHeader();  // to prevent coalescing rightPart with fBlock
                coalescAndPut(rightPart, fBlockEnd - (uintptr_t)rightPart);
            }
            // ... and free left part
            if (newB != fBlock) {
                newB->initHeader(); // to prevent coalescing fBlock with newB
                coalescAndPut(fBlock, (uintptr_t)newB - (uintptr_t)fBlock);
            }

            fBlock = newB;
            MALLOC_ASSERT(isAligned(fBlock, slabSize), ASSERT_TEXT);
        } else {
            if (size_t splitSz = fBlock->sizeTmp - num*size) {
                // split block and return free right part
                FreeBlock *splitB = (FreeBlock*)((uintptr_t)fBlock + num*size);
                splitB->initHeader();
                coalescAndPut(splitB, splitSz);
            }
        }
        bkndSync.signal();
        FreeBlock::markBlocks(fBlock, num, size);
    }

    return fBlock;
}

// try to allocate size Byte block from any of slab-aligned spaces.
// needAlignedRes is true if result must be slab-aligned
FreeBlock *Backend::getFromAlignedSpace(int binIdx, int num, size_t size,
                                        bool needAlignedRes, bool wait, int *binLocked)
{
    FreeBlock *fBlock =
        freeAlignedBins.getBlock(binIdx, &bkndSync, num*size, needAlignedRes,
                                 /*alignedBin=*/true, wait, binLocked);

    if (fBlock) {
        if (fBlock->sizeTmp != num*size) { // i.e., need to split the block
            FreeBlock *newAlgnd;
            size_t newSz;

            if (needAlignedRes) {
                newAlgnd = fBlock;
                fBlock = (FreeBlock*)((uintptr_t)newAlgnd + newAlgnd->sizeTmp
                                      - num*size);
                MALLOC_ASSERT(isAligned(fBlock, slabSize), "Invalid free block");
                fBlock->initHeader();
                newSz = newAlgnd->sizeTmp - num*size;
            } else {
                newAlgnd = (FreeBlock*)((uintptr_t)fBlock + num*size);
                newSz = fBlock->sizeTmp - num*size;
                newAlgnd->initHeader();
            }
            coalescAndPut(newAlgnd, newSz);
        }
        bkndSync.signal();
        MALLOC_ASSERT(!needAlignedRes || isAligned(fBlock, slabSize), ASSERT_TEXT);
        FreeBlock::markBlocks(fBlock, num, size);
    }
    return fBlock;
}

void Backend::correctMaxRequestSize(size_t requestSize)
{
    // Find maximal requested size limited by getMaxBinnedSize()
    if (requestSize < getMaxBinnedSize()) {
        for (size_t oldMaxReq = maxRequestedSize;
             requestSize > oldMaxReq && requestSize < getMaxBinnedSize(); ) {
            size_t val = AtomicCompareExchange((intptr_t&)maxRequestedSize,
                                               requestSize, oldMaxReq);
            if (val == oldMaxReq)
                break;
            oldMaxReq = val;
        }
    }
}

inline size_t Backend::getMaxBinnedSize()
{
    return hugePages.wasObserved && !inUserPool()?
        maxBinned_HugePage : maxBinned_SmallPage;
}

bool Backend::askMemFromOS(size_t blockSize, intptr_t startModifiedCnt,
                           int *lockedBinsThreshold,
                           int numOfLockedBins, bool *largeBinsUpdated)
{
    size_t maxBinSize = 0;

    // Another thread is modifying backend while we can't get the block.
    // Wait while it leaves and re-do the scan
    // before trying other ways to extend the backend.
    if (bkndSync.waitTillSignalled(startModifiedCnt)
        // semaphore is protecting adding more more memory from OS
        || memExtendingSema.wait())
        return true;

    if (startModifiedCnt != bkndSync.getNumOfMods()) {
        memExtendingSema.signal();
        return true;
    }
    // To keep objects below maxBinnedSize, region must be larger then that.
    // So trying to balance between too small regions (that leads to
    // fragmentation) and too large ones (that leads to excessive address
    // space consumption). If region is "quite large", allocate only one,
    // to prevent fragmentation. It supposedly doesn't hurt performance,
    // because the object requested by user is large.
    const size_t maxBinned = getMaxBinnedSize();
    const size_t regSz_sizeBased = blockSize>=maxBinned?
        blockSize : alignUp(4*maxRequestedSize, 1024*1024);
    if (blockSize == slabSize || blockSize == numOfSlabAllocOnMiss*slabSize
        || regSz_sizeBased < maxBinned) {
        for (unsigned idx=0; idx<4; idx++) {
            size_t binSize = addNewRegion(maxBinned, /*exact=*/false);
            if (!binSize)
                break;
            if (binSize > maxBinSize)
                maxBinSize = binSize;
        }
    } else {
        // if huge pages enabled and blockSize>=maxBinned, rest of space up to
        // huge page alignment is unusable, because single user object sits
        // in an region.
        *largeBinsUpdated = true;
        maxBinSize = addNewRegion(regSz_sizeBased, /*exact=*/true);
    }
    memExtendingSema.signal();
    askMemFromOSCounter.OSasked();

    // When blockSize >= maxBinnedSize, and getRawMem failed
    // for this allocation, allocation in bins
    // is our last chance to fulfil the request.
    // Sadly, size is larger then max bin, so have to give up.
    if (maxBinSize && maxBinSize < blockSize)
        return false;

    if (!maxBinSize) { // no regions have been added, try to clean cache
        if (extMemPool->hardCachesCleanup())
            *largeBinsUpdated = true;
        else {
            if (bkndSync.waitTillSignalled(startModifiedCnt))
                return true;
            // OS can't give us more memory, but we have some in locked bins
            if (*lockedBinsThreshold && numOfLockedBins) {
                *lockedBinsThreshold = 0;
                return true;
            }
            return false;
        }
    }
    return true;
}

// try to allocate size Byte block in available bins
// needAlignedRes is true if result must be slab-aligned
FreeBlock *Backend::genericGetBlock(int num, size_t size, bool needAlignedRes)
{
    // after (soft|hard)CachesCleanup we can get memory in large bins,
    // while after addNewRegion only in slab-aligned bins. This flag
    // is for large bins update status.
    bool largeBinsUpdated = true;
    FreeBlock *block = NULL;
    const size_t totalReqSize = num*size;
    const int nativeBin = sizeToBin(totalReqSize);
    // If we found 2 or less locked bins, it's time to ask more memory from OS.
    // But nothing can be asked from fixed pool. And we prefer wait, not ask
    // for more memory, if block is quite large.
    int lockedBinsThreshold = extMemPool->fixedPool || size>=maxBinned_SmallPage? 0 : 2;

    correctMaxRequestSize(totalReqSize);
    scanCoalescQ(/*forceCoalescQDrop=*/false);

    for (;;) {
        const intptr_t startModifiedCnt = bkndSync.getNumOfMods();
        int numOfLockedBins;

        for (;;) {
            numOfLockedBins = 0;

            // TODO: try different bin search order
            if (needAlignedRes) {
                if (!block)
                    for ( int i=freeAlignedBins.getMinNonemptyBin(nativeBin);
                          i<freeBinsNum; i=freeAlignedBins.getMinNonemptyBin(i+1) ){
                        block = getFromAlignedSpace(i, num, size, /*needAlignedRes=*/true, /*wait=*/false, &numOfLockedBins);
                        if (block) break;
                    }
                if (!block && largeBinsUpdated)
                    for ( int i=freeLargeBins.getMinNonemptyBin(nativeBin);
                          i<freeBinsNum; i=freeLargeBins.getMinNonemptyBin(i+1) ){
                        block = getFromBin(i, num, size, /*needAlignedRes=*/true, &numOfLockedBins);
                        if (block) break;
                    }
            } else {
                if (!block && largeBinsUpdated)
                    for ( int i=freeLargeBins.getMinNonemptyBin(nativeBin);
                          i<freeBinsNum; i=freeLargeBins.getMinNonemptyBin(i+1) ){
                        block = getFromBin(i, num, size, /*needAlignedRes=*/false, &numOfLockedBins);
                        if (block) break;
                    }
                if (!block)
                    for ( int i=freeAlignedBins.getMinNonemptyBin(nativeBin);
                          i<freeBinsNum; i=freeAlignedBins.getMinNonemptyBin(i+1) ){
                        block = getFromAlignedSpace(i, num, size, /*needAlignedRes=*/false, /*wait=*/false, &numOfLockedBins);
                        if (block) break;
                    }
            }
            if (block || numOfLockedBins<=lockedBinsThreshold)
                break;
        }
        if (block)
            break;

        largeBinsUpdated = scanCoalescQ(/*forceCoalescQDrop=*/true);
        largeBinsUpdated = extMemPool->softCachesCleanup() || largeBinsUpdated;
        if (!largeBinsUpdated) {
            if (!askMemFromOS(totalReqSize, startModifiedCnt, &lockedBinsThreshold,
                              numOfLockedBins, &largeBinsUpdated))
                return NULL;
        }
    }
    return block;
}

LargeMemoryBlock *Backend::getLargeBlock(size_t size)
{
    LargeMemoryBlock *lmb =
        (LargeMemoryBlock*)genericGetBlock(1, size, /*needAlignedRes=*/false);
    if (lmb) {
        lmb->unalignedSize = size;
        if (extMemPool->mustBeAddedToGlobalLargeBlockList())
            extMemPool->lmbList.add(lmb);
    }
    return lmb;
}

void *Backend::getBackRefSpace(size_t size, bool *rawMemUsed)
{
    // This block is released only at shutdown, so it can prevent
    // a entire region releasing when it's received from the backend,
    // so prefer getRawMemory using.
    if (void *ret = getRawMemory(size, /*hugePages=*/false)) {
        *rawMemUsed = true;
        return ret;
    }
    void *ret = genericGetBlock(1, size, /*needAlignedRes=*/false);
    if (ret) *rawMemUsed = false;
    return ret;
}

void Backend::putBackRefSpace(void *b, size_t size, bool rawMemUsed)
{
    if (rawMemUsed)
        freeRawMemory(b, size);
    // ignore not raw mem, as it released on region releasing
}

void Backend::removeBlockFromBin(FreeBlock *fBlock)
{
    if (fBlock->myBin != Backend::NO_BIN) {
        if (fBlock->aligned)
            freeAlignedBins.lockRemoveBlock(fBlock->myBin, fBlock);
        else
            freeLargeBins.lockRemoveBlock(fBlock->myBin, fBlock);
    }
}

void Backend::genericPutBlock(FreeBlock *fBlock, size_t blockSz)
{
    bkndSync.consume();
    coalescAndPut(fBlock, blockSz);
    bkndSync.signal();
}

void AllLargeBlocksList::add(LargeMemoryBlock *lmb)
{
    MallocMutex::scoped_lock scoped_cs(largeObjLock);
    lmb->gPrev = NULL;
    lmb->gNext = loHead;
    if (lmb->gNext)
        lmb->gNext->gPrev = lmb;
    loHead = lmb;
}

void AllLargeBlocksList::remove(LargeMemoryBlock *lmb)
{
    MallocMutex::scoped_lock scoped_cs(largeObjLock);
    if (loHead == lmb)
        loHead = lmb->gNext;
    if (lmb->gNext)
        lmb->gNext->gPrev = lmb->gPrev;
    if (lmb->gPrev)
        lmb->gPrev->gNext = lmb->gNext;
}

void AllLargeBlocksList::removeAll(Backend *backend)
{
    LargeMemoryBlock *next, *lmb = loHead;
    loHead = NULL;

    for (; lmb; lmb = next) {
        next = lmb->gNext;
        // nothing left to AllLargeBlocksList::remove
        lmb->gNext = lmb->gPrev = NULL;
        backend->returnLargeObject(lmb);
    }
}

void Backend::putLargeBlock(LargeMemoryBlock *lmb)
{
    if (extMemPool->mustBeAddedToGlobalLargeBlockList())
        extMemPool->lmbList.remove(lmb);
    genericPutBlock((FreeBlock *)lmb, lmb->unalignedSize);
}

void Backend::returnLargeObject(LargeMemoryBlock *lmb)
{
    removeBackRef(lmb->backRefIdx);
    putLargeBlock(lmb);
    STAT_increment(getThreadId(), ThreadCommonCounters, freeLargeObj);
}

void Backend::releaseRegion(MemRegion *memRegion)
{
    {
        MallocMutex::scoped_lock lock(regionListLock);
        if (regionList == memRegion)
            regionList = memRegion->next;
        if (memRegion->next)
            memRegion->next->prev = memRegion->prev;
        if (memRegion->prev)
            memRegion->prev->next = memRegion->next;
    }
    freeRawMem(memRegion, memRegion->allocSz);
}

// coalesce fBlock with its neighborhood
FreeBlock *Backend::doCoalesc(FreeBlock *fBlock, MemRegion **mRegion)
{
    FreeBlock *resBlock = fBlock;
    size_t resSize = fBlock->sizeTmp;
    MemRegion *memRegion = NULL;

    fBlock->markCoalescing(resSize);
    resBlock->blockInBin = false;

    // coalesing with left neighbor
    size_t leftSz = fBlock->trySetLeftUsed(GuardedSize::COAL_BLOCK);
    if (leftSz != GuardedSize::LOCKED) {
        if (leftSz == GuardedSize::COAL_BLOCK) {
            coalescQ.putBlock(fBlock);
            return NULL;
        } else {
            FreeBlock *left = fBlock->leftNeig(leftSz);
            size_t lSz = left->trySetMeUsed(GuardedSize::COAL_BLOCK);
            if (lSz <= GuardedSize::MAX_LOCKED_VAL) {
                fBlock->setLeftFree(leftSz); // rollback
                coalescQ.putBlock(fBlock);
                return NULL;
            } else {
                MALLOC_ASSERT(lSz == leftSz, "Invalid header");
                left->blockInBin = true;
                resBlock = left;
                resSize += leftSz;
                resBlock->sizeTmp = resSize;
            }
        }
    }
    // coalesing with right neighbor
    FreeBlock *right = fBlock->rightNeig(fBlock->sizeTmp);
    size_t rightSz = right->trySetMeUsed(GuardedSize::COAL_BLOCK);
    if (rightSz != GuardedSize::LOCKED) {
        // LastFreeBlock is on the right side
        if (GuardedSize::LAST_REGION_BLOCK == rightSz) {
            right->setMeFree(GuardedSize::LAST_REGION_BLOCK);
            memRegion = static_cast<LastFreeBlock*>(right)->memRegion;
        } else if (GuardedSize::COAL_BLOCK == rightSz) {
            if (resBlock->blockInBin) {
                resBlock->blockInBin = false;
                removeBlockFromBin(resBlock);
            }
            coalescQ.putBlock(resBlock);
            return NULL;
        } else {
            size_t rSz = right->rightNeig(rightSz)->
                trySetLeftUsed(GuardedSize::COAL_BLOCK);
            if (rSz <= GuardedSize::MAX_LOCKED_VAL) {
                right->setMeFree(rightSz);  // rollback
                if (resBlock->blockInBin) {
                    resBlock->blockInBin = false;
                    removeBlockFromBin(resBlock);
                }
                coalescQ.putBlock(resBlock);
                return NULL;
            } else {
                MALLOC_ASSERT(rSz == rightSz, "Invalid header");
                removeBlockFromBin(right);
                resSize += rightSz;

                // Is LastFreeBlock on the right side of right?
                FreeBlock *nextRight = right->rightNeig(rightSz);
                size_t nextRightSz = nextRight->
                    trySetMeUsed(GuardedSize::COAL_BLOCK);
                if (nextRightSz > GuardedSize::MAX_LOCKED_VAL) {
                    if (nextRightSz == GuardedSize::LAST_REGION_BLOCK)
                        memRegion = static_cast<LastFreeBlock*>(nextRight)->memRegion;

                    nextRight->setMeFree(nextRightSz);
                }
            }
        }
    }
    if (memRegion) {
        MALLOC_ASSERT((uintptr_t)memRegion + memRegion->allocSz >=
                      (uintptr_t)right + sizeof(LastFreeBlock), ASSERT_TEXT);
        MALLOC_ASSERT((uintptr_t)memRegion < (uintptr_t)resBlock, ASSERT_TEXT);
        *mRegion = memRegion;
    } else
        *mRegion = NULL;
    resBlock->sizeTmp = resSize;
    return resBlock;
}

void Backend::coalescAndPutList(FreeBlock *list, bool forceCoalescQDrop)
{
    FreeBlock *helper;
    MemRegion *memRegion;

    for (;list; list = helper) {
        bool addToTail = false;
        helper = list->nextToFree;
        FreeBlock *toRet = doCoalesc(list, &memRegion);
        if (!toRet)
            continue;

        if (memRegion && memRegion->blockSz == toRet->sizeTmp
            && !extMemPool->fixedPool) {
            if (extMemPool->regionsAreReleaseable()) {
                // release the region, because there is no used blocks in it
                if (toRet->blockInBin)
                    removeBlockFromBin(toRet);
                releaseRegion(memRegion);
                continue;
            } else // add block from empty region to end of bin,
                addToTail = true; // preserving for exact fit
        }
        size_t currSz = toRet->sizeTmp;
        int bin = sizeToBin(currSz);
        bool toAligned = toAlignedBin(toRet, currSz);
        bool needAddToBin = true;

        if (toRet->blockInBin) {
            // is it stay in same bin?
            if (toRet->myBin == bin && toRet->aligned == toAligned)
                needAddToBin = false;
            else {
                toRet->blockInBin = false;
                removeBlockFromBin(toRet);
            }
        }

        // not stay in same bin, or bin-less, add it
        if (needAddToBin) {
            toRet->prev = toRet->next = toRet->nextToFree = NULL;
            toRet->myBin = NO_BIN;

            // If the block is too small to fit in any bin, keep it bin-less.
            // It's not a leak because the block later can be coalesced.
            if (currSz >= minBinnedSize) {
                toRet->sizeTmp = currSz;
                IndexedBins *target = toAligned? &freeAlignedBins : &freeLargeBins;
                if (forceCoalescQDrop) {
                    target->tryAddBlock(bin, toRet, addToTail);
                } else if (!target->tryAddBlock(bin, toRet, addToTail)) {
                    coalescQ.putBlock(toRet);
                    continue;
                }
            }
            toRet->sizeTmp = 0;
        }
        // Free (possibly coalesced) free block.
        // Adding to bin must be done before this point,
        // because after a block is free it can be coalesced, and
        // using its pointer became unsafe.
        // Remember that coalescing is not done under any global lock.
        toRet->setMeFree(currSz);
        toRet->rightNeig(currSz)->setLeftFree(currSz);
    }
}

// Coalesce fBlock and add it back to a bin;
// processing delayed coalescing requests.
void Backend::coalescAndPut(FreeBlock *fBlock, size_t blockSz)
{
    fBlock->sizeTmp = blockSz;
    fBlock->nextToFree = NULL;

    coalescAndPutList(fBlock, /*forceCoalescQDrop=*/false);
}

bool Backend::scanCoalescQ(bool forceCoalescQDrop)
{
    FreeBlock *currCoalescList = coalescQ.getAll();

    if (currCoalescList)
        coalescAndPutList(currCoalescList, forceCoalescQDrop);
    return currCoalescList;
}

FreeBlock *Backend::findBlockInRegion(MemRegion *region)
{
    FreeBlock *fBlock;
    size_t blockSz;
    uintptr_t fBlockEnd,
        lastFreeBlock = (uintptr_t)region + region->allocSz - sizeof(LastFreeBlock);

    if (region->exact) {
        fBlock = (FreeBlock *)alignUp((uintptr_t)region + sizeof(MemRegion),
                                      largeObjectAlignment);
        fBlockEnd = lastFreeBlock;
    } else { // right bound is slab-aligned, keep LastFreeBlock after it
        fBlock = (FreeBlock *)((uintptr_t)region + sizeof(MemRegion));
        fBlockEnd = alignDown(lastFreeBlock, slabSize);
    }
    if (fBlockEnd <= (uintptr_t)fBlock)
        return NULL; // allocSz is too small
    blockSz = fBlockEnd - (uintptr_t)fBlock;
    // TODO: extend getSlabBlock to support degradation, i.e. getting less blocks
    // then requested, and then relax this check
    // (now all or nothing is implemented, check according to this)
    if (blockSz < numOfSlabAllocOnMiss*slabSize)
        return NULL;

    region->blockSz = blockSz;
    return fBlock;
}

// startUseBlock adds free block to a bin, the block can be used and
// even released after this, so the region must be added to regionList already
void Backend::startUseBlock(MemRegion *region, FreeBlock *fBlock)
{
    size_t blockSz = region->blockSz;
    fBlock->initHeader();
    fBlock->setMeFree(blockSz);

    LastFreeBlock *lastBl = static_cast<LastFreeBlock*>(fBlock->rightNeig(blockSz));
    lastBl->initHeader();
    lastBl->setMeFree(GuardedSize::LAST_REGION_BLOCK);
    lastBl->setLeftFree(blockSz);
    lastBl->myBin = NO_BIN;
    lastBl->memRegion = region;

    unsigned targetBin = sizeToBin(blockSz);
    if (!region->exact && toAlignedBin(fBlock, blockSz)) {
        freeAlignedBins.addBlock(targetBin, fBlock, blockSz, /*addToTail=*/false);
    } else {
        freeLargeBins.addBlock(targetBin, fBlock, blockSz, /*addToTail=*/false);
    }
}

size_t Backend::addNewRegion(size_t rawSize, bool exact)
{
    // to guarantee that header is not overwritten in used blocks
    MALLOC_ASSERT(sizeof(BlockMutexes) <= sizeof(BlockI), ASSERT_TEXT);
    // to guarantee that block length is not conflicting with
    // special values of GuardedSize
    MALLOC_ASSERT(FreeBlock::minBlockSize > GuardedSize::MAX_SPEC_VAL, ASSERT_TEXT);
    // "exact" means that not less than rawSize for block inside the region.
    // Reserve space for region header, worst case alignment
    // and last block mark.
    if (exact)
        rawSize += sizeof(MemRegion) + largeObjectAlignment 
                +  FreeBlock::minBlockSize + sizeof(LastFreeBlock);

    MemRegion *region = (MemRegion*)getRawMem(rawSize);
    if (!region) return 0;
    if (rawSize < sizeof(MemRegion)) {
        if (!extMemPool->fixedPool)
            freeRawMem(region, rawSize);
        return 0;
    }

    region->exact = exact;
    region->allocSz = rawSize;
    FreeBlock *fBlock = findBlockInRegion(region);
    if (!fBlock) {
        if (!extMemPool->fixedPool)
            freeRawMem(region, rawSize);
        return 0;
    }
    // adding to global list of all regions
    {
        region->prev = NULL;
        MallocMutex::scoped_lock lock(regionListLock);
        region->next = regionList;
        regionList = region;
        if (regionList->next)
            regionList->next->prev = regionList;
    }
    // copy it here, as just after starting to use region it might be released
    size_t blockSz = region->blockSz;

    startUseBlock(region, fBlock);
    bkndSync.pureSignal();
    return blockSz;
}

void Backend::reset()
{
    MemRegion *curr;

    MALLOC_ASSERT(extMemPool->userPool(), "Only user pool can be reset.");
    // no active threads are allowed in backend while reset() called
    verify();

    freeLargeBins.reset();
    freeAlignedBins.reset();

    for (curr = regionList; curr; curr = curr->next) {
        FreeBlock *fBlock = findBlockInRegion(curr);
        MALLOC_ASSERT(fBlock, "A memory region unexpectedly got smaller");
        startUseBlock(curr, fBlock);
    }
}

bool Backend::destroy()
{
    // no active threads are allowed in backend while destroy() called
    verify();
    while (regionList) {
        MemRegion *helper = regionList->next;
        if (inUserPool())
            (*extMemPool->rawFree)(extMemPool->poolId, regionList,
                                   regionList->allocSz);
        else {
            freeRawMemory(regionList, regionList->allocSz);
        }
        regionList = helper;
    }
    return true;
}

void Backend::IndexedBins::verify()
{
    for (int i=0; i<freeBinsNum; i++) {
        for (FreeBlock *fb = freeBins[i].head; fb; fb=fb->next) {
            uintptr_t mySz = fb->myL.value;
            MALLOC_ASSERT(mySz>GuardedSize::MAX_SPEC_VAL, ASSERT_TEXT);
            FreeBlock *right = (FreeBlock*)((uintptr_t)fb + mySz);
            suppress_unused_warning(right);
            MALLOC_ASSERT(right->myL.value<=GuardedSize::MAX_SPEC_VAL, ASSERT_TEXT);
            MALLOC_ASSERT(right->leftL.value==mySz, ASSERT_TEXT);
            MALLOC_ASSERT(fb->leftL.value<=GuardedSize::MAX_SPEC_VAL, ASSERT_TEXT);
        }
    }
}

// For correct operation, it must be called when no other threads
// is changing backend.
void Backend::verify()
{
#if MALLOC_DEBUG
    scanCoalescQ(/*forceCoalescQDrop=*/false);

    freeLargeBins.verify();
    freeAlignedBins.verify();
#endif // MALLOC_DEBUG
}

#if __TBB_MALLOC_BACKEND_STAT
size_t Backend::Bin::countFreeBlocks()
{
    size_t cnt = 0;
    {
        MallocMutex::scoped_lock lock(tLock);
        for (FreeBlock *fb = head; fb; fb = fb->next)
            cnt++;
    }
    return cnt;
}

void Backend::IndexedBins::reportStat(FILE *f)
{
    size_t totalSize = 0;

    for (int i=0; i<Backend::freeBinsNum; i++)
        if (size_t cnt = freeBins[i].countFreeBlocks()) {
            totalSize += cnt*Backend::binToSize(i);
            fprintf(f, "%d:%lu ", i, cnt);
        }
    fprintf(f, "\ttotal size %lu KB", totalSize/1024);
}

void Backend::reportStat(FILE *f)
{
    int regNum = 0;

    scanCoalescQ(/*forceCoalescQDrop=*/false);

    {
        MallocMutex::scoped_lock lock(regionListLock);
        for (MemRegion *curr = regionList; curr; curr = curr->next)
            regNum++;
    }
    fprintf(f, "%d regions\nlarge ", regNum);
    freeLargeBins.reportStat(f);
    fprintf(f, "\naligned ");
    freeAlignedBins.reportStat(f);
    fprintf(f, "\n");
}
#endif // __TBB_MALLOC_BACKEND_STAT

} } // namespaces
