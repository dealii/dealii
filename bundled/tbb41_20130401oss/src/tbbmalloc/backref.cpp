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
#include <string.h>
#include <new>        /* for placement new */

namespace rml {
namespace internal {


/********* backreferences ***********************/
/* Each slab block and each large memory object header contains BackRefIdx
 * that points out in some BackRefBlock which points back to this block or header.
 */
struct BackRefBlock : public BlockI {
    BackRefBlock *nextForUse;     // the next in the chain of blocks with free items
    FreeObject   *bumpPtr;        // bump pointer moves from the end to the beginning of the block
    FreeObject   *freeList;
    // list of all blocks that were allocated from raw mem (i.e., not from backend)
    BackRefBlock *nextRawMemBlock;
    int           allocatedCount; // the number of objects allocated
    int           myNum;          // the index in the parent array
    MallocMutex   blockMutex;
    bool          addedToForUse;  // this block is already added to the listForUse chain

    BackRefBlock(BackRefBlock *blockToUse, int num) :
        nextForUse(NULL), bumpPtr((FreeObject*)((uintptr_t)blockToUse + slabSize - sizeof(void*))),
        freeList(NULL), nextRawMemBlock(NULL), allocatedCount(0), myNum(num),
        addedToForUse(false) {
        memset(&blockMutex, 0, sizeof(MallocMutex));
        // index in BackRefMaster must fit to uint16_t
        MALLOC_ASSERT(!(myNum >> 16), ASSERT_TEXT);
    }

    // TODO: take into account VirtualAlloc granularity
    static const int bytes = slabSize;
};

// max number of backreference pointers in slab block
static const int BR_MAX_CNT = (BackRefBlock::bytes-sizeof(BackRefBlock))/sizeof(void*);

struct BackRefMaster {
/* A slab block can hold up to ~2K back pointers to slab blocks or large objects,
 * so it can address at least 32MB. The array of 64KB holds 8K pointers
 * to such blocks, addressing ~256 GB.
 */
    static const size_t bytes = 64*1024;
    static const int dataSz;
/* space is reserved for master table and 4 leaves
   taking into account VirtualAlloc allocation granularity */
    static const int leaves = 4;
    static const size_t masterSize = BackRefMaster::bytes+leaves*BackRefBlock::bytes;

    Backend       *backend;
    BackRefBlock  *active;         // if defined, use it for allocations
    BackRefBlock  *listForUse;     // the chain of data blocks with free items
    BackRefBlock  *allRawMemBlocks;
    intptr_t       lastUsed;       // index of the last used block
    bool           rawMemUsed;
    BackRefBlock  *backRefBl[1];   // the real size of the array is dataSz

    BackRefBlock *findFreeBlock();
    void          addBackRefBlockToList(BackRefBlock *bl);
    void          addEmptyBackRefBlock(BackRefBlock *newBl);
};

const int BackRefMaster::dataSz
    = 1+(BackRefMaster::bytes-sizeof(BackRefMaster))/sizeof(BackRefBlock*);

static MallocMutex backRefMutex;
static BackRefMaster *backRefMaster;

bool initBackRefMaster(Backend *backend)
{
    bool rawMemUsed;
    BackRefMaster *master =
        (BackRefMaster*)backend->getBackRefSpace(BackRefMaster::masterSize,
                                                 &rawMemUsed);
    if (! master)
        return false;
    master->backend = backend;
    master->listForUse = master->allRawMemBlocks = NULL;
    master->rawMemUsed = rawMemUsed;
    master->lastUsed = -1;
    for (int i=0; i<BackRefMaster::leaves; i++) {
        BackRefBlock *bl = (BackRefBlock *)((uintptr_t)master + BackRefMaster::bytes + i*BackRefBlock::bytes);
        master->addEmptyBackRefBlock(bl);
        if (i)
            master->addBackRefBlockToList(bl);
        else // active leaf is not needed in listForUse
            master->active = bl;
    }
    // backRefMaster is read in getBackRef, so publish it in consistent state
    FencedStore((intptr_t&)backRefMaster, (intptr_t)master);
    return true;
}

void destroyBackRefMaster(Backend *backend)
{
    if (backRefMaster) { // Is initBackRefMaster() called?
        for (BackRefBlock *curr=backRefMaster->allRawMemBlocks; curr; ) {
            BackRefBlock *next = curr->nextRawMemBlock;
            // allRawMemBlocks list is only for raw mem blocks
            backend->putBackRefSpace(curr, BackRefBlock::bytes, /*rawMemUsed=*/true);
            curr = next;
        }
        backend->putBackRefSpace(backRefMaster, BackRefMaster::masterSize,
                                 backRefMaster->rawMemUsed);
    }
}

void BackRefMaster::addBackRefBlockToList(BackRefBlock *bl)
{
    bl->nextForUse = listForUse;
    listForUse = bl;
    bl->addedToForUse = true;
}

void BackRefMaster::addEmptyBackRefBlock(BackRefBlock *newBl)
{
    intptr_t nextLU = lastUsed+1;
    memset((char*)newBl+sizeof(BackRefBlock), 0,
           BackRefBlock::bytes-sizeof(BackRefBlock));
    new (newBl) BackRefBlock(newBl, nextLU);
    backRefBl[nextLU] = newBl;
    // lastUsed is read in getBackRef, and access to backRefBl[lastUsed]
    // is possible only after checking backref against current lastUsed
    FencedStore(lastUsed, nextLU);
}

BackRefBlock *BackRefMaster::findFreeBlock()
{
    if (active->allocatedCount < BR_MAX_CNT)
        return active;

    if (listForUse) {                                   // use released list
        active = listForUse;
        listForUse = listForUse->nextForUse;
        MALLOC_ASSERT(active->addedToForUse, ASSERT_TEXT);
        active->addedToForUse = false;
    } else if (lastUsed-1 < backRefMaster->dataSz) {    // allocate new data node
        bool rawMemUsed;
        BackRefBlock *newBl =
            (BackRefBlock*)backend->getBackRefSpace(BackRefBlock::bytes, &rawMemUsed);
        if (!newBl) return NULL;
        backRefMaster->addEmptyBackRefBlock(newBl);
        if (rawMemUsed) {
            newBl->nextRawMemBlock = backRefMaster->allRawMemBlocks;
            backRefMaster->allRawMemBlocks = newBl;
        } else
            newBl->nextRawMemBlock = NULL;
        active = newBl;
    } else  // no free blocks, give up
        return NULL;
    return active;
}

void *getBackRef(BackRefIdx backRefIdx)
{
    // !backRefMaster means no initialization done, so it can't be valid memory
    // see addEmptyBackRefBlock for fences around lastUsed
    if (!FencedLoad((intptr_t&)backRefMaster)
        || backRefIdx.getMaster() > FencedLoad(backRefMaster->lastUsed)
        || backRefIdx.getOffset() >= BR_MAX_CNT)
        return NULL;
    return *(void**)((uintptr_t)backRefMaster->backRefBl[backRefIdx.getMaster()]
                     + sizeof(BackRefBlock)+backRefIdx.getOffset()*sizeof(void*));
}

void setBackRef(BackRefIdx backRefIdx, void *newPtr)
{
    MALLOC_ASSERT(backRefIdx.getMaster()<=backRefMaster->lastUsed && backRefIdx.getOffset()<BR_MAX_CNT,
                  ASSERT_TEXT);
    *(void**)((uintptr_t)backRefMaster->backRefBl[backRefIdx.getMaster()]
              + sizeof(BackRefBlock) + backRefIdx.getOffset()*sizeof(void*)) = newPtr;
}

BackRefIdx BackRefIdx::newBackRef(bool largeObj)
{
    BackRefBlock *blockToUse;
    void **toUse;
    BackRefIdx res;

    do {
        { // global lock taken to find a block
            MallocMutex::scoped_lock lock(backRefMutex);

            MALLOC_ASSERT(backRefMaster, ASSERT_TEXT);
            blockToUse = backRefMaster->findFreeBlock();
            if (!blockToUse)
                return BackRefIdx();
        }
        toUse = NULL;
        { // the block is locked to find a reference
            MallocMutex::scoped_lock lock(blockToUse->blockMutex);

            if (blockToUse->freeList) {
                toUse = (void**)blockToUse->freeList;
                blockToUse->freeList = blockToUse->freeList->next;
                MALLOC_ASSERT(!blockToUse->freeList ||
                              ((uintptr_t)blockToUse->freeList>=(uintptr_t)blockToUse
                               && (uintptr_t)blockToUse->freeList <
                               (uintptr_t)blockToUse + slabSize), ASSERT_TEXT);
            } else if (blockToUse->allocatedCount < BR_MAX_CNT) {
                toUse = (void**)blockToUse->bumpPtr;
                blockToUse->bumpPtr =
                    (FreeObject*)((uintptr_t)blockToUse->bumpPtr - sizeof(void*));
                if (blockToUse->allocatedCount == BR_MAX_CNT-1) {
                    MALLOC_ASSERT((uintptr_t)blockToUse->bumpPtr
                                  < (uintptr_t)blockToUse+sizeof(BackRefBlock),
                                  ASSERT_TEXT);
                    blockToUse->bumpPtr = NULL;
                }
            }
            if (toUse)
                blockToUse->allocatedCount++;
        } // end of lock scope
    } while (!toUse);
    res.master = blockToUse->myNum;
    uintptr_t offset =
        ((uintptr_t)toUse - ((uintptr_t)blockToUse + sizeof(BackRefBlock)))/sizeof(void*);
    // Is offset too big?
    MALLOC_ASSERT(!(offset >> 15), ASSERT_TEXT);
    res.offset = offset;
    if (largeObj) res.largeObj = largeObj;

    return res;
}

void removeBackRef(BackRefIdx backRefIdx)
{
    MALLOC_ASSERT(!backRefIdx.isInvalid(), ASSERT_TEXT);
    MALLOC_ASSERT(backRefIdx.getMaster()<=backRefMaster->lastUsed
                  && backRefIdx.getOffset()<BR_MAX_CNT, ASSERT_TEXT);
    BackRefBlock *currBlock = backRefMaster->backRefBl[backRefIdx.getMaster()];
    FreeObject *freeObj = (FreeObject*)((uintptr_t)currBlock + sizeof(BackRefBlock)
                                        + backRefIdx.getOffset()*sizeof(void*));
    MALLOC_ASSERT(((uintptr_t)freeObj>(uintptr_t)currBlock &&
                   (uintptr_t)freeObj<(uintptr_t)currBlock + slabSize), ASSERT_TEXT);
    {
        MallocMutex::scoped_lock lock(currBlock->blockMutex);

        freeObj->next = currBlock->freeList;
        MALLOC_ASSERT(!freeObj->next ||
                      ((uintptr_t)freeObj->next > (uintptr_t)currBlock
                       && (uintptr_t)freeObj->next <
                       (uintptr_t)currBlock + slabSize), ASSERT_TEXT);
        currBlock->freeList = freeObj;
        currBlock->allocatedCount--;
    }
    // TODO: do we need double-check here?
    if (!currBlock->addedToForUse && currBlock!=backRefMaster->active) {
        MallocMutex::scoped_lock lock(backRefMutex);

        if (!currBlock->addedToForUse && currBlock!=backRefMaster->active)
            backRefMaster->addBackRefBlockToList(currBlock);
    }
}

/********* End of backreferences ***********************/

} // namespace internal
} // namespace rml

